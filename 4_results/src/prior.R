library(argparse)
library(Matrix)
library(dplyr)
library(GpGp)
library(Rcpp)
library(fastsparse)

source(Sys.getenv('UTILS_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--basis-vectors')
parser$add_argument('--perturbations')
parser$add_argument('--output')
args <- parser$parse_args()

basis_vectors <- fst::read_fst(args$basis_vectors)
perturbations <- fst::read_fst(args$perturbations)
alpha_mean <- rep(0, nrow(basis_vectors))

with(basis_vectors, {
  is_ocean_inventory <<- inventory == 'ocean'
  is_small_bio_region <<- (
    region %in% c(sprintf('Region%02d', 12 : 22), 'RegionNZ')
    & startsWith(as.character(inventory), 'bio')
  )
  is_climatology <<- component != 'residual'
})

alpha_precision <- diag(ifelse(
  is_ocean_inventory,
  # Ocean prior is actually set by CAMS-based prior below
  1 / 100,
  1
))

ocean_perturbations <- perturbations %>%
  filter(
    inventory == 'ocean',
    component == 'residual'
  ) %>%
  add_basis_vector(basis_vectors)

ocean_locations <- ocean_perturbations %>%
  select(longitude, latitude, time) %>%
  distinct() %>%
  mutate(location_index = 1 : n())

ocean_location_matrix <- as.matrix(
  ocean_locations %>%
    mutate(time_i = as.integer(factor(time))) %>%
    select(longitude, latitude, time_i)
)

ocean_location_nn <- find_ordered_nn(ocean_location_matrix, m = 48)

# CAMS prior is 0.1 gC/m^2/d on 8 day grid cells with e-folding time 1 month.
# Grid cell size is roughly similar. To convert to monthly, we imagine
# averaging 4 grid cells (32 days) as follows:
Q <- ou_precision(
  c(1, 1, 1) / 31 * 8,
  1,
  rep(1 / 0.1 ^ 2, 4),
  rep(1 / 0.1 ^ 2, 3)
)
X <- rbind(c(1, 1, 1, 1))
cams_equivalent_sd <- sqrt((X %*% solve(Q, t(X)))[1])

G_C_M2_D_TO_KG_M2_S <- 44.01 / 12.01 / 1e3 / 24 / 60 / 60
cams_equivalent_variance <- (G_C_M2_D_TO_KG_M2_S * cams_equivalent_sd) ^ 2
ocean_Sigma_L_inv <- vecchia_Linv(
  c(
    0.999 * cams_equivalent_variance,
    # e-folding length 1000km
    1000 / 6371,
    # e-folding time 1 month
    1,
    # Nugget
    0.001 * cams_equivalent_variance
  ),
  'exponential_spheretime',
  ocean_location_matrix,
  ocean_location_nn
)

X_ocean_constraint <- with(
  ocean_perturbations %>%
    left_join(ocean_locations, by = c('longitude', 'latitude', 'time')),
  sparseMatrix(
  i = location_index,
  j = as.integer(basis_vector),
  x = ocean_perturbations$value,
  dims = c(max(location_index), nlevels(basis_vector))
))

L_inv_X_ocean_constraint_parts <- parallel::mclapply(
  which(colSums(abs(X_ocean_constraint)) > 0),
  function(j) {
    part <- Linv_mult(
      ocean_Sigma_L_inv,
      X_ocean_constraint[, j],
      ocean_location_nn
    )
    non_zero <- which(part != 0)
    list(
      i = non_zero,
      j = rep(j, length(non_zero)),
      x = part[non_zero]
    )
  },
  mc.cores = get_cores()
)
L_inv_X_ocean_constraint <- sparseMatrix(
  i = do.call(c, lapply(L_inv_X_ocean_constraint_parts, getElement, 'i')),
  j = do.call(c, lapply(L_inv_X_ocean_constraint_parts, getElement, 'j')),
  x = do.call(c, lapply(L_inv_X_ocean_constraint_parts, getElement, 'x')),
  dims = c(nrow(X_ocean_constraint), ncol(X_ocean_constraint))
)

ocean_constraint_Xt_Q_X <- as.matrix(crossprod(L_inv_X_ocean_constraint))
alpha_precision_final <- alpha_precision + as.matrix(ocean_constraint_Xt_Q_X)

# Exclude some regions
alpha_to_exclude <- (
  (is_small_bio_region & is_climatology)
  | (is_ocean_inventory & is_climatology)
)
diag(alpha_precision_final) <- ifelse(
  alpha_to_exclude,
  Inf,
  diag(alpha_precision_final)
)

saveRDS(list(
  mean = alpha_mean,
  precision = alpha_precision_final
), args$output)
