library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(fastsparse)
library(WoodburyMatrix)

source(Sys.getenv('UTILS_PARTIAL'))

nlm_best_of <- function(fn, k, max_attempts, best_of, ...) {
  good_attempts <- 0
  best_result <- list(minimum = Inf)
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch({
      nlm(fn, rnorm(k), ...)
    }, error = function(e) {
      print(e)
      NULL
    })

    if (is.null(result) || result$code %in% (3 : 5)) next

    good_attempts <- good_attempts + 1
    if (result$minimum < best_result$minimum) {
      best_result <- result
    }
    if (good_attempts == best_of) break
  }
  if (attempt == max_attempts) stop('Max attempts exceeded')

  best_result
}

log_posterior <- function(
  residual,
  time_since,
  error,
  rho,
  ell,
  gamma
) {
  log_likelihood <- 0
  for (j in seq_along(rho)) {
    for (i in seq_along(residual[[j]])) {
      diff_time_since <- diff(time_since[[j]][[i]])
      # HACK(mgnb): some IS sites have repeated times; separate them minimally
      diff_time_since[diff_time_since == 0] <- 1 / 24
      stopifnot(all(diff_time_since != 0))
      precisions <- 1 / error[[j]][[i]] ^ 2
      cross_precisions <- head(sqrt(precisions), -1) * tail(sqrt(precisions), -1)

      if (rho[j] == 1) {
        Q <- ou_precision(
          diff_time_since,
          1 / ell[j],
          precisions * gamma,
          cross_precisions * gamma
        )
        log_likelihood <- log_likelihood + (
          0.5 * determinant(Q)$modulus
          - 0.5 * as.vector(crossprod(residual[[j]][[i]], Q %*% residual[[j]][[i]]))
        )
      } else {
        A <- FastDiagonal(
          x = (gamma / (1 - rho[j])) * precisions
        )
        B <- ou_precision(
          diff_time_since,
          1 / ell[j],
          precisions * gamma / rho[j],
          cross_precisions * gamma / rho[j]
        )
        O <- TridiagonalMatrix(
          A@x + B@major,
          B@minor
        )
        Sigma <- WoodburyMatrix(
          A,
          B,
          O = O,
          symmetric = TRUE
        )
        log_likelihood <- log_likelihood + dwnorm(residual[[j]][[i]], covariance = Sigma, log = TRUE)
      }
    }
  }

  (
    log_likelihood
    + dgamma(gamma, shape = GAMMA_PRIOR$shape, rate = GAMMA_PRIOR$rate, log = TRUE)
    + sum(dgamma(ell, shape = ELL_PRIOR$shape, rate = ELL_PRIOR$rate, log = TRUE))
  )
}

parser <- ArgumentParser()
parser$add_argument('--residuals')
parser$add_argument('--observations')
parser$add_argument('--output')
args <- parser$parse_args()

log_debug('Loading residuals')
residuals <- fst::read_fst(args$residuals) %>%
  left_join(
    fst::read_fst(args$observations) %>%
      select(
        observation_id,
        overall_observation_mode,
        observation_group,
        hyperparameter_group,
        time,
        error
      ),
    by = 'observation_id'
  )

log_debug('Performing fits')
# NOTE(mgnb): each hyperparameter_group has its own rho, ell, and is nested
# within a fit_group, which has its own gamma
fits <- residuals %>%
  arrange(time) %>%
  mutate(
    fit_group = if_else(
      overall_observation_mode == 'IS',
      '3_IS',
      as.character(observation_group)
    )
  ) %>%
  group_by(fit_group) %>%
  group_modify(~ {
    log_trace('Fitting {.y$fit_group}')
    is_oco2 <- .y$fit_group %in% c('1_LNLG', '2_OG')

    unit <- if (is_oco2) 'mins' else 'days'

    # NOTE(mgnb): split first by hyperparameter group, then obs group
    group_parts <- .x %>%
      mutate(
        time_since = as.double(time - min(time), unit = unit),
      ) %>%
      group_by(hyperparameter_group) %>%
      group_map(~ {
        .x %>%
          group_by(observation_group) %>%
          group_map(~ .x)
      })

    n_groups <- length(group_parts)
    residual_parts <- lapply(group_parts, lapply, getElement, 'residual')
    time_since_parts <- lapply(group_parts, lapply, getElement, 'time_since')
    error_parts <- lapply(group_parts, lapply, getElement, 'error')

    fit <- if (is_oco2) {
      nlm_best_of(function(params) {
        gamma <- exp(params[1])
        rho <- 1 / (1 + exp(-params[2 : (n_groups + 1)]))
        ell <- exp(tail(params, n_groups))

        -log_posterior(
          residual_parts,
          time_since_parts,
          error_parts,
          gamma = gamma,
          rho = rho,
          ell = ell
        )
      }, 1 + 2 * n_groups, max_attempts = 10, best_of = 3, print.level = 0)
    } else {
      nlm_best_of(function(params) {
        gamma <- exp(params[1])
        rho <- rep(1, n_groups)
        ell <- exp(tail(params, n_groups))

        -log_posterior(
          residual_parts,
          time_since_parts,
          error_parts,
          gamma = gamma,
          rho = rho,
          ell = ell
        )
      }, 1 + n_groups, max_attempts = 10, best_of = 3, print.level = 0)
    }

    if (is_oco2) {
      gamma <- exp(fit$estimate[1])
      rho <- 1 / (1 + exp(-fit$estimate[2 : (n_groups + 1)]))
      ell <- exp(tail(fit$estimate, n_groups))
    } else {
      gamma <- exp(fit$estimate[1])
      rho <- rep(1, n_groups)
      ell <- exp(tail(fit$estimate, n_groups))
    }

    data.frame(
      n = nrow(.x),
      hyperparameter_group = sort(unique(.x$hyperparameter_group)),
      rho = rho,
      ell = ell,
      gamma = gamma,
      ell_unit = unit
    )
  }) %>%
  ungroup() %>%
  select(hyperparameter_group, rho, ell, gamma, ell_unit)

fits <- bind_rows(
  fits,
  fits %>%
    filter(hyperparameter_group == 'surface') %>%
    mutate(hyperparameter_group = 'lauder')
)

log_debug('Saving to {args$output}')
fst::write_fst(fits, args$output)

log_debug('Done')
