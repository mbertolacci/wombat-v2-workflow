library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(fastsparse)
library(WoodburyMatrix)
library(Matrix)
library(Rcpp)
library(fst)

rcpp_cache_dir <- Sys.getenv('RCPP_CACHE_DIR')
options(rcpp.cache.dir = if (rcpp_cache_dir == '') tempdir() else rcpp_cache_dir)

source(Sys.getenv('UTILS_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))
sourceCpp(Sys.getenv('HMC_EXACT_CPP_PARTIAL'))

slice <- function(w, min_w = 0, max_evaluations = 100) {
  learning <- TRUE
  sum <- 0
  n <- 0

  function(x0, log_f, learn = FALSE, include_n_evaluations = FALSE) {
    n_evaluations <- 1
    log_f_count <- function(x) {
      output <- log_f(x)
      n_evaluations <<- n_evaluations + 1
      stopifnot(n_evaluations <= max_evaluations)
      output
    }

    if (learn) {
      if (!learning) learning <<- TRUE
    } else {
      if (learning) {
        if (n > 0) w <<- max(min_w, sum / n)
        learning <<- FALSE
      }
    }

    y <- log(runif(1)) + log_f_count(x0)
    if (y == -Inf) stop('x0 value outside the support')
    U <- runif(1)
    L <- x0 - w * U
    R <- L + w

    while (y < log_f_count(L)) L <- L - w
    while (y < log_f_count(R)) R <- R + w

    repeat {
      x1 <- runif(1, L, R)
      if (y < log_f_count(x1)) break
      if (x1 < x0) L <- x1 else R <- x1
    }

    if (learning) {
      sum <<- sum + abs(x1 - x0) # as suggested, sec 4.4
      n <<- n + 1
    }

    if (include_n_evaluations) {
      list(n_evaluations = n_evaluations, w = w, sample = x1)
    } else {
      x1
    }
  }
}

dmvnorm <- function(x, mean, covariance, precision, log = FALSE) {
  if (missing(mean)) mean <- 0

  z <- x - mean
  output <- if (!missing(covariance)) {
    as.numeric(
      - 0.5 * (length(x) * log(2 * pi))
      - 0.5 * determinant(covariance, logarithm = TRUE)$modulus
      - 0.5 * crossprod(z, solve(covariance, z))
    )
  } else {
    if (is.matrix(precision)) {
      chol_precision <- chol(precision)
      as.numeric(
        - 0.5 * (length(x) * log(2 * pi))
        + sum(log(diag(chol_precision)))
        - 0.5 * crossprod(chol_precision %*% z)
      )
    } else {
      as.numeric(
        - 0.5 * (length(x) * log(2 * pi))
        + 0.5 * determinant(precision, logarithm = TRUE)$modulus
        - 0.5 * crossprod(z, precision %*% z)
      )
    }
  }

  if (log) output else exp(output)
}

parser <- ArgumentParser()
parser$add_argument('--overall-observation-mode', nargs = '+')
parser$add_argument('--control', nargs = '+')
parser$add_argument('--constraints')
parser$add_argument('--component-name', nargs = '+')
parser$add_argument('--component-parts', nargs = '+')
parser$add_argument('--component-transport-matrix', nargs = '+')
parser$add_argument('--observations')
parser$add_argument('--hyperparameter-estimates')
parser$add_argument('--basis-vectors')
parser$add_argument('--prior')
parser$add_argument('--output')
args <- parser$parse_args()

log_debug('Loading control')
control <- bind_rows(lapply(args$control, read_fst)) %>%
  mutate(observation_id = droplevels(observation_id))

log_debug('Loading observations')
observations <- read_fst(args$observations) %>%
  filter(
    assimilate %in% c(1, 3),
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
stopifnot(!anyNA(observations$error))

basis_vectors <- read_fst(args$basis_vectors)
prior <- readRDS(args$prior)

n_all_alpha <- nrow(basis_vectors)
# NOTE(mgnb): values of alpha fixed to zero are given infinite precision
alpha_to_include <- is.finite(diag(prior$precision)) & with(
  basis_vectors,
  !(inventory == 'bio_resp_tot' & component %in% c('intercept', 'trend'))
)

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
H_component_parts <- lapply(part_indices, function(part_i) {
  log_debug('Loading {args$component_transport_matrix[part_i]}')
  name_i <- args$component_name[part_i]
  observations_i <- observations %>%
    filter(component_name == name_i)
  n_observations_i <- observations_i %>%
    nrow() %>%
    as.double()

  n_i <- n_observations_i * n_all_alpha
  log_trace('Total size to load is {n_i} = {n_observations_i} x {n_all_alpha}')

  fn <- pipe(sprintf('lz4 -v %s -', args$component_transport_matrix[part_i]), 'rb')
  H_vec <- readBin(fn, 'double', n_i)
  close(fn)
  output <- matrix(H_vec, nrow = n_observations_i)
  gc()
  output[, alpha_to_include]
})
gc()

hyperparameter_groups <- observations %>%
  group_by(hyperparameter_group) %>%
  group_modify(~ tibble(
    component_names = list(unique(.x$component_name))
  ))
hyperparameter_group_indices <- seq_len(nrow(hyperparameter_groups))

log_debug('Setting up MCMC')
H_parts <- lapply(hyperparameter_group_indices, function(i) {
  component_names_i <- hyperparameter_groups$component_names[[i]]
  hyperparameter_group_i <- hyperparameter_groups$hyperparameter_group[i]

  H_parts_i <- lapply(component_names_i, function(component_name_i_j) {
    observations_component_i_j <- observations %>%
      filter(
        component_name == component_name_i_j
      )
    component_index <- which(args$component_name == component_name_i_j)
    H_component_parts[[component_index]][
      observations_component_i_j$hyperparameter_group == hyperparameter_group_i,
    ]
  })
  if (length(H_parts_i) == 1) {
    H_parts_i[[1]]
  } else {
    do.call(rbind, H_parts_i)
  }
})
rm(H_component_parts)

observation_parts <- lapply(hyperparameter_group_indices, function(i) {
  component_names_i <- hyperparameter_groups$component_names[[i]]
  hyperparameter_group_i <- hyperparameter_groups$hyperparameter_group[i]

  lapply(component_names_i, function(component_name_i_j) {
    observations %>%
      filter(
        component_name == component_name_i_j,
        hyperparameter_group == hyperparameter_group_i
      )
  }) %>%
    bind_rows()
})

offset_parts <- lapply(observation_parts, getElement, 'offset')

hyperparameter_estimates <- read_fst(args$hyperparameter_estimates)

Sigma_epsilon_parts <- lapply(hyperparameter_group_indices, function(i) {
  observations_i <- observation_parts[[i]]

  is_precision <- NA
  parts <- observations_i %>%
    group_by(observation_group, hyperparameter_group) %>%
    group_map(~ {
      parameters <- .y %>%
        left_join(hyperparameter_estimates, by = 'hyperparameter_group')

      precisions <- 1 / .x$error ^ 2
      cross_precisions <- head(sqrt(precisions), -1) * tail(sqrt(precisions), -1)
      diff_time <- diff(as.double(
        .x$time - .x$time[1],
        unit = parameters$ell_unit
      ))
      diff_time[diff_time == 0] <- 1 / 24

      if (parameters$rho == 1) {
        is_precision <<- TRUE
        ou_precision(
          diff_time,
          1 / parameters$ell,
          precisions,
          cross_precisions
        )
      } else {
        is_precision <<- FALSE
        A <- FastDiagonal(
          x = (1 / (1 - parameters$rho)) * precisions
        )
        B <- ou_precision(
          diff_time,
          1 / parameters$ell,
          precisions / parameters$rho,
          cross_precisions / parameters$rho
        )
        O <- TridiagonalMatrix(
          A@x + B@major,
          B@minor
        )
        WoodburyMatrix(
          A,
          B,
          O = O,
          symmetric = TRUE
        )
      }
    })
  if (!is_precision) {
    # NOTE(mgnb): implicitly, this is a WoodburyMatrix
    output <- parts[[1]]
    attr(output, 'is_Q') <- FALSE
    output
  } else {
    # NOTE(mgnb): implicitly, this is a list of TridiagonalMatrix's
    output <- bdiag_tridiagonal(parts)
    attr(output, 'is_Q') <- TRUE
    output
  }
})

Ht_Q_epsilon_H_parts <- lapply(hyperparameter_group_indices, function(i) {
  log_trace('Computing Ht_Q_epsilon_H for {hyperparameter_groups$hyperparameter_group[i]} (n = {nrow(H_parts[[i]])})')
  gc()
  if (is(Sigma_epsilon_parts[[i]], 'SWoodburyMatrix')) {
    # HACK(mgnb): this is a slower but more memory efficient way of calculating
    # the solve, necessary because the matrices can get so huge
    a <- Sigma_epsilon_parts[[i]]
    b <- H_parts[[i]]
    rhs <- a@A %*% b
    rhs <- solve(a@O, rhs)
    gc()
    rhs <- (-a@A) %*% rhs
    gc()
    rhs <- rhs + a@A %*% b
    gc()

    as.matrix(crossprod(H_parts[[i]], rhs))
  } else {
    if (attr(Sigma_epsilon_parts[[i]], 'is_Q')) {
      as.matrix(crossprod(H_parts[[i]], Sigma_epsilon_parts[[i]] %*% H_parts[[i]]))
    } else {
      as.matrix(crossprod(H_parts[[i]], solve(Sigma_epsilon_parts[[i]], H_parts[[i]])))
    }
  }
})

Ht_Q_epsilon_y_parts <- lapply(hyperparameter_group_indices, function(i) {
  log_trace('Computing Ht_Q_epsilon_y for {hyperparameter_groups$hyperparameter_group[i]}')
  if (attr(Sigma_epsilon_parts[[i]], 'is_Q')) {
    as.vector(
      crossprod(H_parts[[i]], Sigma_epsilon_parts[[i]] %*% offset_parts[[i]])
    )
  } else {
    as.vector(
      crossprod(H_parts[[i]], solve(Sigma_epsilon_parts[[i]], offset_parts[[i]]))
    )
  }
})
rm(H_parts)
gc()

yt_Q_epsilon_y_parts <- lapply(hyperparameter_group_indices, function(i) {
  log_trace('Computing yt_Q_epsilon_y for {hyperparameter_groups$hyperparameter_group[i]}')
  if (attr(Sigma_epsilon_parts[[i]], 'is_Q')) {
    as.vector(
      crossprod(offset_parts[[i]], Sigma_epsilon_parts[[i]] %*% offset_parts[[i]])
    )
  } else {
    as.vector(
      crossprod(offset_parts[[i]], solve(Sigma_epsilon_parts[[i]], offset_parts[[i]]))
    )
  }
})

alpha_prior_mean <- prior$mean[alpha_to_include]
alpha_prior_precision_base <- prior$precision[alpha_to_include, alpha_to_include]

bio_regions <- sort(unique(basis_vectors$region))[1 : 11]
bio_inventories <- c('bio_assim', 'bio_resp_tot')

n_alpha <- length(alpha_prior_mean)
n_bio_regions <- length(bio_regions)

climatology_bio_assim_linear_indices <- with(
  basis_vectors[alpha_to_include, ],
  which(inventory == 'bio_assim' & component %in% c('intercept', 'trend'))
)
climatology_bio_assim_linear_diagonals <- cbind(
  climatology_bio_assim_linear_indices,
  climatology_bio_assim_linear_indices
)
climatology_bio_non_linear_indices <- with(basis_vectors[alpha_to_include, ], list(
  bio_assim = which(inventory == 'bio_assim' & !(component %in% c('intercept', 'trend', 'residual'))),
  bio_resp_tot = which(inventory == 'bio_resp_tot' & !(component %in% c('intercept', 'trend', 'residual')))
))
climatology_bio_non_linear_diagonals <- cbind(
  # NOTE(mgnb): the rbind cause these to alternate between assim/resp
  as.vector(do.call(rbind, climatology_bio_non_linear_indices)),
  as.vector(do.call(rbind, climatology_bio_non_linear_indices))
)
climatology_bio_non_linear_off_diagonals <- rbind(
  cbind(
    climatology_bio_non_linear_indices$bio_assim,
    climatology_bio_non_linear_indices$bio_resp_tot
  ),
  cbind(
    climatology_bio_non_linear_indices$bio_resp_tot,
    climatology_bio_non_linear_indices$bio_assim
  )
)

residual_bio_indices_i <- lapply(bio_regions, function(region_i) {
  which(with(basis_vectors[alpha_to_include, ], {
    region == region_i & inventory != 'ocean' & component == 'residual'
  }))
})

n_times <- length(residual_bio_indices_i[[1]]) / 2

get_alpha_prior_precision <- function(
  rho_bio_clim, w_bio_clim,
  kappa_bio_resid, rho_bio_resid, w_bio_resid
) {
  output <- alpha_prior_precision_base

  # Bio climatology process
  output[climatology_bio_assim_linear_diagonals] <- w_bio_clim[1]
  Q_pair <- solve(rbind(
    c(1 / w_bio_clim[1], rho_bio_clim / sqrt(prod(w_bio_clim))),
    c(rho_bio_clim / sqrt(prod(w_bio_clim)), 1 / w_bio_clim[2])
  ))
  output[climatology_bio_non_linear_diagonals] <- diag(Q_pair)
  output[climatology_bio_non_linear_off_diagonals] <- Q_pair[1, 2]

  # Residual process
  Q_ar <- ar1_Q(n_times, kappa_bio_resid, sparse = FALSE)
  Q_cross <- solve(rbind(
    c(1 / w_bio_resid[1], rho_bio_resid / sqrt(prod(w_bio_resid))),
    c(rho_bio_resid / sqrt(prod(w_bio_resid)), 1 / w_bio_resid[2])
  ))
  Q <- kronecker(Q_ar, Q_cross)
  for (i in seq_along(bio_regions)) {
    indices <- residual_bio_indices_i[[i]]
    output[indices, indices] <- Q
  }

  output
}

log_pdf_climatology_bio <- function(rho, w) {
  if (rho < -1 || rho > 1) return(-Inf)
  if (any(w < 0)) return(-Inf)
  Q_pair <- solve(rbind(
    c(1 / w[1], rho / sqrt(prod(w))),
    c(rho / sqrt(prod(w)), 1 / w[2])
  ))
  indices <- as.vector(do.call(rbind, climatology_bio_non_linear_indices))
  Q <- kronecker(
    diag(length(indices) / 2),
    Q_pair
  )
  sum(dnorm(
    alpha_current[climatology_bio_assim_linear_indices],
    alpha_prior_mean[climatology_bio_assim_linear_indices],
    sd = 1 / sqrt(w[1]),
    log = TRUE
  )) + dmvnorm(
    alpha_current[indices],
    alpha_prior_mean[indices],
    precision = Q,
    log = TRUE
  ) + sum(dgamma(
    w,
    shape = W_PRIOR$shape,
    rate = W_PRIOR$rate,
    log = TRUE
  ))
}

log_pdf_residual_bio <- function(kappa, rho, w, alpha) {
  if (kappa < 0 | kappa > 1) return(-Inf)
  if (rho < -1 | rho > 1) return(-Inf)
  if (any(w < 0)) return(-Inf)

  Q_ar <- ar1_Q(n_times, kappa, sparse = FALSE)
  Q_cross <- solve(rbind(
    c(1 / w[1], rho / sqrt(prod(w))),
    c(rho / sqrt(prod(w)), 1 / w[2])
  ))
  Q <- kronecker(Q_ar, Q_cross)
  log_likelihoods <- sapply(seq_along(bio_regions), function(i) {
    dmvnorm(
      alpha_current[residual_bio_indices_i[[i]]],
      alpha_prior_mean[residual_bio_indices_i[[i]]],
      precision = Q,
      log = TRUE
    )
  })

  sum(log_likelihoods) + sum(dgamma(
    w,
    shape = W_PRIOR$shape,
    rate = W_PRIOR$rate,
    log = TRUE
  ))
}

if (!is.null(args$constraints)) {
  constraints <- readRDS(args$constraints)
  F_constraint <- rbind(
    constraints$F_sign,
    constraints$F_residual
  )[, alpha_to_include]
  g_constraint <- c(constraints$g_sign, constraints$g_residual)
}

n_hyperparameter_groups <- nrow(hyperparameter_groups)

alpha_samples <- matrix(NA, nrow = N_MCMC_SAMPLES, ncol = n_alpha)
rho_bio_clim_samples <- rep(NA, N_MCMC_SAMPLES)
w_bio_clim_samples <- matrix(NA, nrow = N_MCMC_SAMPLES, ncol = 2)
kappa_bio_resid_samples <- rep(NA, N_MCMC_SAMPLES)
rho_bio_resid_samples <- rep(NA, N_MCMC_SAMPLES)
w_bio_resid_samples <- matrix(NA, nrow = N_MCMC_SAMPLES, ncol = 2)
gamma_samples <- matrix(NA, nrow = N_MCMC_SAMPLES, ncol = n_hyperparameter_groups)

# NOTE(mgnb): this starting value satisfies the initial condition
alpha_current <- rep(0, n_alpha)
rho_bio_clim_current <- 0
w_bio_clim_current <- c(1, 1)
kappa_bio_resid_current <- 0
rho_bio_resid_current <- 0
w_bio_resid_current <- c(1, 1)
gamma_current <- rep(1,  n_hyperparameter_groups)

alpha_samples[1, ] <- alpha_current
rho_bio_clim_samples[1] <- rho_bio_clim_current
w_bio_clim_samples[1, ] <- w_bio_clim_current
kappa_bio_resid_samples[1] <- kappa_bio_resid_current
rho_bio_resid_samples[1] <- rho_bio_resid_current
w_bio_resid_samples[1, ] <- w_bio_resid_current
gamma_samples[1, ] <- gamma_current

rho_bio_clim_slice <- slice(0.1, max_evaluations = 1000)
w_bio_clim_slice <- list(
  slice(0.1, max_evaluations = 1000),
  slice(0.1, max_evaluations = 1000)
)
kappa_bio_resid_slice <- slice(0.1, max_evaluations = 1000)
rho_bio_resid_slice <- slice(0.1, max_evaluations = 1000)
w_bio_resid_slice <- list(
  slice(0.1, max_evaluations = 1000),
  slice(0.1, max_evaluations = 1000)
)

log_debug('Starting MCMC')
for (iteration in 2 : N_MCMC_SAMPLES) {
  log_trace('[{iteration} / {N_MCMC_SAMPLES}] Sampling alpha')
  alpha_prior_precision_current <- get_alpha_prior_precision(
    rho_bio_clim_current,
    w_bio_clim_current,
    kappa_bio_resid_current,
    rho_bio_resid_current,
    w_bio_resid_current
  )
  alpha_precision <- (
    Reduce(`+`, lapply(seq_len(n_hyperparameter_groups), function(i) {
      gamma_current[i] * Ht_Q_epsilon_H_parts[[i]]
    }))
    + alpha_prior_precision_current
  )
  chol_alpha_precision <- chol(alpha_precision)
  alpha_mean <- chol_solve(
    chol_alpha_precision,
    Reduce(`+`, lapply(seq_len(n_hyperparameter_groups), function(i) {
      gamma_current[i] * Ht_Q_epsilon_y_parts[[i]]
    }))
    + alpha_prior_precision_current %*% alpha_prior_mean
  )
  if (!is.null(args$constraints)) {
    alpha_current <- sampleHmcConstrained(
      rep(0, n_alpha),
      alpha_mean,
      chol_alpha_precision,
      F_constraint,
      g_constraint,
      pi / 2,
      nSamples = 1,
      debug = TRUE,
      bounceLimit = 100000
    )
  } else {
    alpha_current <- as.vector(
      alpha_mean + backsolve(chol_alpha_precision, rnorm(n_alpha))
    )
  }

  log_trace('[{iteration} / {N_MCMC_SAMPLES}] Sampling rho_bio_clim (current = {format(rho_bio_clim_current)})')
  rho_bio_clim_current <- rho_bio_clim_slice(rho_bio_clim_current, function(rho) {
    log_pdf_climatology_bio(rho, w_bio_clim_current)
  }, learn = iteration <= N_MCMC_WARM_UP)

  log_trace('[{iteration} / {N_MCMC_SAMPLES}] Sampling w_bio_clim (current = {paste0(format(w_bio_clim_current), collapse = ", ")})')
  for (i in c(1, 2)) {
    w_bio_clim_current[i] <- w_bio_clim_slice[[i]](w_bio_clim_current[i], function(w_i) {
      w <- w_bio_clim_current
      w[i] <- w_i
      log_pdf_climatology_bio(rho_bio_clim_current, w)
    }, learn = iteration <= N_MCMC_WARM_UP)
  }

  log_trace('[{iteration} / {N_MCMC_SAMPLES}] Sampling kappa_bio_resid (current = {format(kappa_bio_resid_current)})')
  kappa_bio_resid_current <- kappa_bio_resid_slice(kappa_bio_resid_current, function(kappa) {
    log_pdf_residual_bio(kappa, rho_bio_resid_current, w_bio_resid_current, alpha_current)
  }, learn = iteration <= N_MCMC_WARM_UP)

  log_trace('[{iteration} / {N_MCMC_SAMPLES}] Sampling rho_bio_resid (current = {format(rho_bio_resid_current)})')
  rho_bio_resid_current <- rho_bio_resid_slice(rho_bio_resid_current, function(rho) {
    log_pdf_residual_bio(kappa_bio_resid_current, rho, w_bio_resid_current, alpha_current)
  }, learn = iteration <= N_MCMC_WARM_UP)

  log_trace('[{iteration} / {N_MCMC_SAMPLES}] Sampling w_bio_resid (current = {paste0(format(w_bio_resid_current), collapse = ", ")})')
  for (i in c(1, 2)) {
    w_bio_resid_current[i] <- w_bio_resid_slice[[i]](w_bio_resid_current[i], function(w_i) {
      w <- w_bio_resid_current
      w[i] <- w_i
      log_pdf_residual_bio(kappa_bio_resid_current, rho_bio_resid_current, w, alpha_current)
    }, learn = iteration <= N_MCMC_WARM_UP)
  }

  log_trace('[{iteration} / {N_MCMC_SAMPLES}] Sampling gamma (current = {paste0(format(gamma_current), collapse = ", ")})')
  for (i in seq_len(n_hyperparameter_groups)) {
    gamma_current[i] <- rgamma(
      1,
      shape = GAMMA_PRIOR$shape + length(offset_parts[[i]]) / 2,
      rate = GAMMA_PRIOR$rate + (
        yt_Q_epsilon_y_parts[[i]]
        - 2 * sum(alpha_current * Ht_Q_epsilon_y_parts[[i]])
        + crossprod(alpha_current, Ht_Q_epsilon_H_parts[[i]] %*% alpha_current)[1]
      ) / 2
    )
  }

  alpha_samples[iteration, ] <- alpha_current
  rho_bio_clim_samples[iteration] <- rho_bio_clim_current
  w_bio_clim_samples[iteration, ] <- w_bio_clim_current
  kappa_bio_resid_samples[iteration] <- kappa_bio_resid_current
  rho_bio_resid_samples[iteration] <- rho_bio_resid_current
  w_bio_resid_samples[iteration, ] <- w_bio_resid_current
  gamma_samples[iteration, ] <- gamma_current
}

colnames(alpha_samples) <- basis_vectors[alpha_to_include, ]$basis_vector_str
colnames(w_bio_clim_samples) <- bio_inventories
colnames(w_bio_resid_samples) <- bio_inventories
colnames(gamma_samples) <- hyperparameter_groups$hyperparameter_group

log_debug('Saving to {args$output}')
alpha_df <- cbind(
  basis_vectors[alpha_to_include, ],
  data.frame(
    value = colMeans(tail(alpha_samples, N_MCMC_SAMPLES - N_MCMC_WARM_UP))
  )
)
alpha_df$value_samples <- t(tail(alpha_samples, N_MCMC_SAMPLES - N_MCMC_WARM_UP))

output <- list(
  alpha = alpha_samples,
  alpha_df = alpha_df,
  rho_bio_clim = rho_bio_clim_samples,
  w_bio_clim = w_bio_clim_samples,
  kappa_bio_resid = kappa_bio_resid_samples,
  rho_bio_resid = rho_bio_resid_samples,
  w_bio_resid = w_bio_resid_samples,
  gamma = gamma_samples
)

saveRDS(output, args$output)
log_debug('Done')
