library(argparse)
library(dplyr)
library(Rcpp)
library(Matrix)
library(logger)
library(parallel)

source(Sys.getenv('UTILS_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--basis-vectors')
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations')
parser$add_argument('--output')
args <- parser$parse_args()

basis_vectors <- fst::read_fst(args$basis_vectors)
control_emissions <- fst::read_fst(args$control_emissions) %>%
  filter(inventory %in% c('bio_assim', 'bio_resp_tot'))
perturbations <- fst::read_fst(args$perturbations) %>%
  filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
  add_basis_vector(basis_vectors)

# Construct constraint for the cell total
perturbations_region <- perturbations %>%
  mutate(
    value = ifelse(inventory == 'bio_assim', -value, value)
  )
baseline_region_cell <- perturbations_region %>%
  group_by(region, inventory, longitude, latitude, time) %>%
  summarise(value = sum(value)) %>%
  ungroup()
F_cell <- with(
  perturbations_region %>%
    left_join(
      baseline_region_cell %>%
        select(-value) %>%
        mutate(cell_index = 1 : n()),
      by = c('inventory', 'longitude', 'latitude', 'time')
    ),
  sparseMatrix(
    i = cell_index,
    j = as.integer(basis_vector),
    x = value,
    dims = c(nrow(baseline_region_cell), nrow(basis_vectors))
  )
)
g_cell <- baseline_region_cell$value

# Construct constraint for the cell climatology
perturbations_region_clim <- perturbations_region %>%
  filter(component != 'residual')
baseline_region_clim <- perturbations_region %>%
  group_by(region, inventory, longitude, latitude, time) %>%
  summarise(value = sum(value)) %>%
  ungroup()
F_clim <- with(
  perturbations_region %>%
    left_join(
      baseline_region_clim %>%
        select(-value) %>%
        mutate(cell_index = 1 : n()),
      by = c('inventory', 'longitude', 'latitude', 'time')
    ),
  sparseMatrix(
    i = cell_index,
    j = as.integer(basis_vector),
    x = value,
    dims = c(nrow(baseline_region_clim), nrow(basis_vectors))
  )
)
g_clim <- baseline_region_clim$value

# Construct constraint for residual
with(
  basis_vectors %>%
    filter(
      inventory %in% c('bio_assim', 'bio_resp_tot'),
      component == 'residual'
    ),
  {
    baseline_residual <<- data.frame(
      region = region,
      inventory = inventory,
      longitude = NA,
      latitude = NA,
      time = NA,
      value = NA
    )
    F_residual <<- sparseMatrix(
      i = seq_along(basis_vector),
      j = as.integer(basis_vector),
      dims = c(length(basis_vector), nlevels(basis_vector))
    )
    g_residual <<- rep(1, length(basis_vector))
  }
)

F <- rbind(F_cell, F_clim)
baseline <- rbind(baseline_region_cell, baseline_region_clim)
g <- pmax(1e-10, c(g_cell, g_clim))

get_redundant_indices <- function(row_indices, column_indices) {
  F_alt <- slam::as.simple_triplet_matrix(-1e7 * F[row_indices, column_indices])
  g_alt <- 1e7 * g[row_indices]
  n_constraints <- nrow(F_alt)
  n_variables <- ncol(F_alt)

  removed_indices <- NULL
  failed_indices <- NULL
  F_alt_current <- F_alt
  for (index in seq_len(n_constraints)) {
    # cat('\r', index, '/', n_constraints, '/', nrow(F_alt_current), '/', length(removed_indices))
    g_prime <- g_alt
    g_prime[index] <- g_alt[index] + 0.1

    if (length(removed_indices) > 0) {
      g_prime_current <- g_prime[-removed_indices]
    } else {
      g_prime_current <- g_prime
    }
    output <- Rglpk::Rglpk_solve_LP(
      obj = F_alt[index, ],
      mat = F_alt_current,
      dir = rep('<', nrow(F_alt) - length(removed_indices)),
      rhs = g_prime_current,
      bounds = list(
        lower = list(ind = 1 : n_variables, val = rep(-Inf, n_variables)),
        upper = list(ind = 1 : n_variables, val = rep(Inf, n_variables))
      ),
      max = TRUE,
      control = list(verbose = FALSE, presolve = FALSE, canonicalize_status = TRUE)
    )
    if (output$status != 0) {
      failed_indices <- c(failed_indices, index)
      next
    }
    if (output$optimum <= g_alt[index]) {
      removed_indices <- c(removed_indices, index)
      F_alt_current <- F_alt[-removed_indices, ]
    }
  }

  if (length(removed_indices) > 0) {
    list(
      removed_indices = row_indices[removed_indices],
      failed_indices = row_indices[failed_indices]
    )
  } else {
    list(
      removed_indices = NULL,
      failed_indices = row_indices[failed_indices]
    )
  }
}

log_info('Starting first round of elimination')
# Eliminate constraints in each grid cell
parts_round_1 <- perturbations %>%
  distinct(inventory, longitude, latitude, region)
redundant_indices_parts_round_1 <- pbmcapply::pbmclapply(seq_len(nrow(parts_round_1)), function(part_index) {
  row_indices <- baseline %>%
    mutate(index = 1 : n()) %>%
    filter(
      inventory == parts_round_1$inventory[part_index],
      longitude == parts_round_1$longitude[part_index],
      latitude == parts_round_1$latitude[part_index],
      region == parts_round_1$region[part_index]
    ) %>%
    pull(index)
  column_indices <- basis_vectors %>%
    mutate(index = 1 : n()) %>%
    filter(
      inventory == parts_round_1$inventory[part_index],
      region == parts_round_1$region[part_index]
    ) %>%
    pull(index)

  get_redundant_indices(row_indices, column_indices)
}, mc.cores = get_cores(), ignore.interactive = TRUE)

redundant_indices_round_1 <- do.call(c, lapply(redundant_indices_parts_round_1, getElement, 'removed_indices'))
failed_indices_round_1 <- do.call(c, lapply(redundant_indices_parts_round_1, getElement, 'failed_indices'))

log_info('First round eliminated {length(redundant_indices_round_1)} rows')

log_info('Starting second round of elimination')
baseline_round_2 <- baseline %>%
  mutate(index = 1 : n()) %>%
  group_by(inventory, region) %>%
  mutate(
    cluster = if (n() <= 1000) 1 else {
      round((1 : n()) / 1000)
    }
  )

parts_round_2 <- baseline_round_2 %>%
  distinct(inventory, region, cluster)
redundant_indices_parts_round_2 <- pbmcapply::pbmclapply(seq_len(nrow(parts_round_2)), function(part_index) {
  row_indices <- baseline_round_2 %>%
    filter(
      inventory == parts_round_2$inventory[part_index],
      region == parts_round_2$region[part_index],
      cluster == parts_round_2$cluster[part_index]
    ) %>%
    pull(index)
  column_indices <- basis_vectors %>%
    mutate(index = 1 : n()) %>%
    filter(
      inventory == parts_round_2$inventory[part_index],
      region == parts_round_2$region[part_index]
    ) %>%
    pull(index)

  # Remove any indices considered redundant in round 1
  row_indices <- setdiff(row_indices, redundant_indices_round_1)

  get_redundant_indices(row_indices, column_indices)
}, mc.cores = get_cores(), ignore.interactive = TRUE)

redundant_indices_round_2 <- do.call(c, lapply(redundant_indices_parts_round_2, getElement, 'removed_indices'))
failed_indices_round_2 <- do.call(c, lapply(redundant_indices_parts_round_2, getElement, 'failed_indices'))

log_info('Second round eliminated {length(redundant_indices_round_2)} rows')

redundant_indices <- sort(unique(c(redundant_indices_round_1, redundant_indices_round_2)))

log_info('Saving to {args$output}')
output <- list(
  F_sign_full = F,
  g_sign_full = g,
  baseline_sign_full = baseline,
  redundant_indices_sign = redundant_indices,
  F_sign = F[-redundant_indices, ],
  g_sign = g[-redundant_indices],
  baseline_sign = baseline[-redundant_indices, ],
  F_residual = F_residual,
  g_residual = g_residual
)
saveRDS(output, args$output)
log_info('Done')
