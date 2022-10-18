source(Sys.getenv('UTILS_PARTIAL'))

library(argparse)
library(dplyr, warn.conflicts = FALSE)

parser <- ArgumentParser()
parser$add_argument('--samples')
parser$add_argument('--output')
args <- parser$parse_args()

samples <- readRDS(args$samples)

print_row <- function(name1, samples1, name2, samples2) {
  samples1 <- samples1[(N_MCMC_WARM_UP + 1) : N_MCMC_SAMPLES]
  samples2 <- samples2[(N_MCMC_WARM_UP + 1) : N_MCMC_SAMPLES]
  cat(sprintf(
    '%s & %.3f & %.3f & %.3f & %s & %.3f & %.3f & %.3f \\\\\n',
    name1,
    quantile(samples1, probs = 0.025),
    quantile(samples1, probs = 0.500),
    quantile(samples1, probs = 0.975),
    name2,
    quantile(samples2, probs = 0.025),
    quantile(samples2, probs = 0.500),
    quantile(samples2, probs = 0.975)
  ))
}

sink(args$output)
cat(
  '\\begin{tabular}{llll||llll}
  \\hline\\
  Variable & 2.5\\% & 50\\% & 97.5\\% & Variable & 2.5\\% & 50\\% & 97.5\\% \\\\ \\hline
')
print_row(
  '$\\tau_\\text{gpp}^\\beta$', samples$w_bio_clim[, 1],
  '$\\kappa_\\text{bio}^\\epsilon$', samples$kappa_bio_resid
)
print_row(
  '$\\tau_\\text{resp}^\\beta$', samples$w_bio_clim[, 2],
  '$(\\gamma_\\text{OCO-2 land}^Z)^{-1}$', 1 / samples$gamma[, 1]
)
print_row(
  '$\\tau_\\text{gpp}^\\epsilon$', samples$w_bio_resid[, 1],
  '$(\\gamma_\\text{Aircraft}^Z)^{-1}$', 1 / samples$gamma[, 2]
)
print_row(
  '$\\tau_\\text{resp}^\\epsilon$', samples$w_bio_resid[, 2],
  '$(\\gamma_\\text{Surface}^Z)^{-1}$', 1 / samples$gamma[, 3]
)
print_row(
  '$\\rho_\\text{gpp,resp}^\\beta$', samples$rho_bio_clim,
  '$(\\gamma_\\text{Tower}^Z)^{-1}$', 1 / samples$gamma[, 4]
)
print_row(
  '$\\rho_\\text{gpp,resp}^\\epsilon$', samples$rho_bio_resid,
  '$(\\gamma_\\text{Shipboard}^Z)^{-1}$', 1 / samples$gamma[, 5]
)
cat('\\hline\n\\end{tabular}')
sink(NULL)
