library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(patchwork)
library(lubridate, warn.conflicts = FALSE)
library(tidyr)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--climatology-by-region')
parser$add_argument('--regions')
parser$add_argument('--sib4-climatology-assim')
parser$add_argument('--output')
args <- parser$parse_args()

climatology_by_region <- readRDS(args$climatology_by_region)

if (args$regions == 'global-zonal') {
  climatology_by_region <- climatology_by_region %>%
    filter(region == 'Global' | grepl('°', as.character(region)))
} else if (args$regions == 'land-transcoms') {
  climatology_by_region <- climatology_by_region %>%
    filter(region %in% sprintf('T%02d', 1 : 11))
} else {
  stop('regions argument not supported')
}

with_nc_file(list(fn = args$sib4_climatology_assim), {
  origin <- ymd_hms(ncdf4::ncatt_get(fn, 0, 'origin')$value)
})

times <- seq(as.Date('2015-01-01'), as.Date('2020-12-01'), by = '1 month')
X_climatology <- climatology_design_matrix(
  time = times,
  origin = origin,
  time_alignment = 'start',
  time_resolution = 'monthly',
  harmonics = 3
)

trend_start <- min(X_climatology[, 2])
trend_end <- max(X_climatology[, 2])

climatology_by_region_harmonic <- climatology_by_region %>%
  filter(!(variable %in% c('intercept', 'trend'))) %>%
  separate(
    variable,
    c('part', 'harmonic', 'coefficient'),
    sep = '_'
  ) %>%
  group_by(output, harmonic, inventory, region, part) %>%
  summarise(
    start = t(value_samples[1, ] + trend_start * value_samples[2, ]),
    end = t(value_samples[1, ] + trend_end * value_samples[2, ])
  ) %>%
  group_by(output, harmonic, inventory, region) %>%
  summarise(
    amplitude_start = t(sqrt(start[1, ] ^ 2 + start[2, ] ^ 2)),
    amplitude_end = t(sqrt(end[1, ] ^ 2 + end[2, ] ^ 2)),
    phase_start = t(atan(-start[2, ] / start[1, ])),
    phase_end = t(atan(-end[2, ] / end[1, ])),
    .groups = 'drop'
  ) %>%
  mutate(
    amplitude_diff = amplitude_end - amplitude_start,
    amplitude_percent = ifelse(
      amplitude_start == 0,
      0,
      100 * (amplitude_end - amplitude_start) / amplitude_start
    ),
    phase_raw_diff = phase_end - phase_start,
    phase_diff = (phase_raw_diff + pi / 2) %% pi - pi / 2,
    phase_diff_days = -ifelse(
      amplitude_start == 0,
      0,
      365.25 * phase_diff / (2 * pi)
    )
  ) %>%
  mutate(
    phase_diff_days_q025 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.025),
    phase_diff_days_q500 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.500),
    phase_diff_days_q975 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.975),
    amplitude_diff_q025 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.025),
    amplitude_diff_q500 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.500),
    amplitude_diff_q975 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.975),
    amplitude_percent_q025 = matrixStats::rowQuantiles(amplitude_percent, probs = 0.025),
    amplitude_percent_q500 = matrixStats::rowQuantiles(amplitude_percent, probs = 0.500),
    amplitude_percent_q975 = matrixStats::rowQuantiles(amplitude_percent, probs = 0.975),
    amplitude_percent_year_q025 = 100 * sign(amplitude_percent_q025) * (((1 + abs(amplitude_percent_q025) / 100) ^ (1 / 6)) - 1),
    amplitude_percent_year_q500 = 100 * sign(amplitude_percent_q500) * (((1 + abs(amplitude_percent_q500) / 100) ^ (1 / 6)) - 1),
    amplitude_percent_year_q975 = 100 * sign(amplitude_percent_q975) * (((1 + abs(amplitude_percent_q975) / 100) ^ (1 / 6)) - 1)
  ) %>%
  mutate(
    region = droplevels(region),
    inventory = c(
      'bio_assim' = 'bioassim',
      'bio_resp_tot' = 'bioresp',
      'nee' = 'NEE'
    )[inventory],
    output = factor(c(
      'Bottom up' = 'Bottom up',
      'Posterior' = 'Posterior'
    )[output])
  ) %>%
  select(
    harmonic,
    output,
    inventory,
    region,
    phase_diff_days_q025,
    phase_diff_days_q975,
    amplitude_diff_q025,
    amplitude_diff_q975,
    amplitude_percent_q025,
    amplitude_percent_q500,
    amplitude_percent_q975,
    amplitude_percent_year_q025,
    amplitude_percent_year_q500,
    amplitude_percent_year_q975
  ) %>%
  arrange(
    harmonic, output, inventory, region
  )

sink(args$output)
knitr::kable(climatology_by_region_harmonic, format = 'simple', digits = 2)
sink(NULL)
