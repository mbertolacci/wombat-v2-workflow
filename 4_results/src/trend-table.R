library(argparse)
library(dplyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--climatology-by-region')
parser$add_argument('--output')
args <- parser$parse_args()

climatology_by_region <- readRDS(args$climatology_by_region)

sink(args$output)
knitr::kable(
  bind_rows(
    climatology_by_region %>%
      filter(variable == 'trend') %>%
      mutate(
        value_q025 = matrixStats::rowQuantiles(value_samples, probs = 0.025),
        value_q975 = matrixStats::rowQuantiles(value_samples, probs = 0.975)
      ) %>%
      select(-value_samples),
  climatology_by_region %>%
    filter(
      variable == 'trend',
      region %in% c('0°-23°N', '23°S-0°')
    ) %>%
    group_by(variable, inventory, output) %>%
    summarise(
      value = sum(value),
      value_q025 = quantile(colSums(value_samples), probs = 0.025),
      value_q975 = quantile(colSums(value_samples), probs = 0.975)
    ) %>%
    mutate(region = 'Tropics')
  ),
  format = 'simple',
  digits = 5
)
sink(NULL)
