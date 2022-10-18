source(Sys.getenv('UTILS_PARTIAL'))

library(argparse)
library(dplyr, warn.conflicts = FALSE)

parser <- ArgumentParser()
parser$add_argument('--hyperparameter-estimates')
parser$add_argument('--output')
args <- parser$parse_args()

hyperparameter_estimates <- fst::read_fst(args$hyperparameter_estimates)

sink(args$output)
cat(sprintf(
  '\\begin{tabular}{lllll}
  \\hline
  OCO-2 land
    & Aircraft
    & Shipboard
    & Surface
    & Tower \\\\
  %.1f seconds
    & %.1f minutes
    & %.1f hours
    & %.1f hours
    & %.1f hours \\\\
  \\hline
\\end{tabular}',
  60 * hyperparameter_estimates$ell[1],
  24 * 60 * hyperparameter_estimates$ell[2],
  24 * hyperparameter_estimates$ell[3],
  24 * hyperparameter_estimates$ell[4],
  24 * hyperparameter_estimates$ell[5]
))
sink(NULL)
