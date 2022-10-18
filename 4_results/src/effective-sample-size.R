source(Sys.getenv('UTILS_PARTIAL'))

library(argparse)

parser <- ArgumentParser()
parser$add_argument('--samples')
parser$add_argument('--output')
args <- parser$parse_args()

samples <- readRDS(args$samples)

ess <- function(x) {
  as.matrix(x)[(N_MCMC_WARM_UP + 1) : N_MCMC_SAMPLES, ] |>
    coda::mcmc() |>
    coda::effectiveSize()
}

sink(args$output)
cat('alpha:\n')
print(quantile(ess(samples$alpha), probs = c(0.05, 0.5)))
for (name in c(
  'w_bio_clim',
  'rho_bio_clim',
  'w_bio_resid',
  'rho_bio_resid',
  'kappa_bio_resid',
  'gamma'
)) {
  cat(name, '\n')
  print(ess(samples[[name]]))
}
sink(NULL)
