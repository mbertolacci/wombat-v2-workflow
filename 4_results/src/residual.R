library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(fastsparse)
library(Matrix)
library(Rcpp)

source(Sys.getenv('UTILS_PARTIAL'))
sourceCpp('scripts/partials/hmc-exact.cpp')

parser <- ArgumentParser()
parser$add_argument('--overall-observation-mode', nargs = '+')
parser$add_argument('--control', nargs = '+')
parser$add_argument('--constraints')
parser$add_argument('--component-name', nargs = '+')
parser$add_argument('--component-parts', nargs = '+')
parser$add_argument('--component-transport-matrix', nargs = '+')
parser$add_argument('--observations')
parser$add_argument('--basis-vectors')
parser$add_argument('--prior')
parser$add_argument('--output')
args <- parser$parse_args()

log_debug('Loading control')
control <- bind_rows(lapply(args$control, fst::read_fst)) %>%
  mutate(observation_id = droplevels(observation_id))

log_debug('Loading observations')
observations <- fst::read_fst(args$observations) %>%
  filter(
    assimilate == 1,
    overall_observation_mode %in% args$overall_observation_mode,
    observation_id %in% control$observation_id
  ) %>%
  arrange(observation_group, time) %>%
  left_join(
    control %>%
      select(observation_id, value_control = value),
    by = 'observation_id'
  ) %>%
  mutate(
    observation_id = factor(
      as.character(observation_id),
      as.character(observation_id)
    ),
    offset = value - value_control
  )

control <- control %>%
  filter(observation_id %in% observations$observation_id) %>%
  mutate(
    observation_id = factor(
      as.character(observation_id),
      levels(observations$observation_id)
    )
  )

stopifnot(all(
  levels(observations$observation_id) == levels(control$observation_id)
))
stopifnot(nlevels(observations$observation_id) == nrow(observations))

basis_vectors <- fst::read_fst(args$basis_vectors)
prior <- readRDS(args$prior)

n_alpha <- nrow(basis_vectors)
# NOTE(mgnb): values of alpha fixed to zero are given infinite precision
alpha_to_include <- is.finite(diag(prior$precision))

part_indices <- seq_along(args$component_name)
observations$component_name <- ''
for (part_i in part_indices) {
  observations$component_name <- ifelse(
    observations$overall_observation_mode %in% strsplit(
      args$component_parts[part_i],
      '|',
      fixed = TRUE
    )[[1]],
    args$component_name[part_i],
    observations$component_name
  )
}

log_debug('Loading transport matrices')
H_parts <- lapply(part_indices, function(part_i) {
  log_debug('Loading {args$component_transport_matrix[part_i]}')
  name_i <- args$component_name[part_i]
  observations_i <- observations %>%
    filter(component_name == name_i)
  n_observations_i <- observations_i %>%
    nrow() %>%
    as.double()

  n_i <- n_observations_i * n_alpha
  log_trace('Total size to load is {n_i} = {n_observations_i} x {n_alpha}')

  fn <- pipe(sprintf('lz4 -v %s -', args$component_transport_matrix[part_i]), 'rb')
  H_vec <- readBin(fn, 'double', n_i)
  close(fn)
  output <- matrix(H_vec, nrow = n_observations_i)
  gc()
  output[, alpha_to_include]
})
gc()

constraints <- readRDS(args$constraints)
F_constraint <- rbind(
  constraints$F_sign,
  constraints$F_residual
)[, alpha_to_include]
g_constraint <- c(constraints$g_sign, constraints$g_residual)

log_debug('Performing inversion')
offset_parts <- lapply(part_indices, function(i) {
  observations %>%
    filter(
      component_name == args$component_name[i]
    ) %>%
    pull(offset)
})
Sigma_epsilon_parts <- lapply(part_indices, function(i) {
  log_trace('Construction Sigma_epsilon for {args$component_name[part_i]}')
  name_i <- args$component_name[i]
  observations_i <- observations %>%
    filter(component_name == name_i)
  FastDiagonal(x = observations_i$error ^ 2)
})
Ht_Q_epsilon_H <- Reduce(`+`, lapply(part_indices, function(part_i) {
  log_trace('Computing Ht_Q_epsilon_H for {args$component_name[part_i]}')
  as.matrix(crossprod(
    H_parts[[part_i]],
    solve(Sigma_epsilon_parts[[part_i]], H_parts[[part_i]])
  ))
}))
Ht_Q_epsilon_y <- Reduce(`+`, lapply(part_indices, function(part_i) {
  log_trace('Computing Ht_Q_epsilon_y for {args$component_name[part_i]}')
  as.vector(crossprod(
    H_parts[[part_i]],
    solve(Sigma_epsilon_parts[[part_i]], offset_parts[[part_i]])
  ))
}))

log_trace('Calculating alpha hat')
alpha_prior_mean <- prior$mean[alpha_to_include]
alpha_prior_precision <- prior$precision[alpha_to_include, alpha_to_include]
alpha_precision <- Ht_Q_epsilon_H + alpha_prior_precision
chol_alpha_precision <- chol(alpha_precision)
alpha_mean <- as.vector(chol_solve(
  chol_alpha_precision,
  Ht_Q_epsilon_y + alpha_prior_precision %*% alpha_prior_mean
))

n_samples <- 100
alpha_samples <- sampleHmcConstrained(
  rep(0, n_alpha),
  alpha_mean,
  chol_alpha_precision,
  F_constraint,
  g_constraint,
  pi / 2,
  nSamples = n_samples,
  debug = TRUE,
  bounceLimit = 5000
)

alpha_hat <- colMeans(alpha_samples)

log_trace('Calculating residuals')
residuals <- bind_rows(lapply(part_indices, function(i) {
  observations %>%
    filter(
      component_name == args$component_name[i]
    ) %>%
    mutate(
      residual = offset - as.vector(H_parts[[i]] %*% alpha_hat)
    ) %>%
    select(observation_id, residual)
}))

log_debug('Saving to {args$output}')
fst::write_fst(residuals, args$output)

log_debug('Done')
