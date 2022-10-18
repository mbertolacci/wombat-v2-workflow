library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(Matrix)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--control-emissions')
parser$add_argument('--region-grid')
parser$add_argument('--samples')
parser$add_argument('--sib4-climatology-assim')
parser$add_argument('--sib4-climatology-resp-tot')
parser$add_argument('--output')
args <- parser$parse_args()

log_debug('Loading inputs')
control_emissions <- fst::read_fst(args$control_emissions)
samples <- readRDS(args$samples)
region_grid <- readRDS(args$region_grid)

fit_region_grid_df <- NULL
for (latitude_index in seq_along(attr(region_grid$Region01, 'latitude'))) {
  for (longitude_index in seq_along(attr(region_grid$Region01, 'longitude'))) {
    found_any <- FALSE
    for (region_name in names(region_grid)) {
      if (region_grid[[region_name]][longitude_index, latitude_index] > 0) {
        fit_region_grid_df <- rbind(fit_region_grid_df, data.frame(
          longitude = attr(region_grid$Region01, 'longitude')[longitude_index],
          latitude = attr(region_grid$Region01, 'latitude')[latitude_index],
          region = region_name
        ))
        found_any <- TRUE
      }
    }
    if (!found_any) {
      fit_region_grid_df <- rbind(fit_region_grid_df, data.frame(
        longitude = attr(region_grid$Region01, 'longitude')[longitude_index],
        latitude = attr(region_grid$Region01, 'latitude')[latitude_index],
        region = 'Region00'
      ))
    }
  }
}

X_fit_region_to_2x25 <- with(
  fit_region_grid_df %>%
    mutate(region = factor(region)),
  sparseMatrix(
    i = seq_along(longitude),
    j = as.integer(region),
    x = rep(1, length(longitude)),
    dims = c(length(longitude), nlevels(region))
  )
)[, 2 : 12]

origin <- ymd_hms(attr(read_gridded_data(
  args$sib4_climatology_assim,
  'coefficient',
  include_time = FALSE,
  include_variable = TRUE
), 'origin'))

trend_prior <- bind_rows(
  read_climatology(args$sib4_climatology_assim) %>%
    mutate(
      value = -value,
      inventory = 'bio_assim'
    ),
  read_climatology(args$sib4_climatology_resp_tot) %>%
    mutate(inventory = 'bio_resp_tot')
) %>%
  filter(variable == 'trend')

log_debug('Computing posterior functionals')
alpha_modified <- samples$alpha_df %>%
  select(
    inventory,
    region,
    variable = component,
    alpha_mean = value,
    alpha_samples = value_samples
  ) %>%
  mutate(
    variable = factor(as.character(variable), levels(trend_prior$variable)),
    alpha_samples = alpha_samples[, floor(seq(1, ncol(alpha_samples), length.out = 200))]
  ) %>%
  filter(variable == 'trend')

trend_posterior <- trend_prior %>%
  group_by(inventory) %>%
  group_modify(~ {
    alpha_samples <- alpha_modified %>%
      filter(inventory == .y$inventory) %>%
      pull(alpha_samples)

    if (nrow(alpha_samples) == 0) {
      alpha_samples <- matrix(0, nrow = 11, ncol = ncol(alpha_samples))
    }

    X_i <- Diagonal(x = .x$value) %*% X_fit_region_to_2x25

    .x %>%
      select(-value) %>%
      mutate(
        value_samples = .x$value + as.matrix(X_i %*% alpha_samples),
        value = rowMeans(value_samples)
      )
  })

output <- bind_rows(
  trend_prior %>%
    mutate(
      output = 'Bottom-up',
      value_samples = matrix(
        value,
        nrow = length(value),
        ncol = ncol(trend_posterior$value_samples)
      )
    ),
  trend_posterior %>% mutate(output = 'Posterior')
) %>%
  group_by(longitude, latitude, output) %>%
  summarise(
    value = sum(value),
    value_samples = t(colSums(value_samples)),
    .groups = 'drop'
  ) %>%
  mutate(
    value_mean = PER_SECONDS_TO_PER_YEAR * ifelse(value == 0, NA, value),
    value_sd = PER_SECONDS_TO_PER_YEAR * ifelse(value == 0, NA, matrixStats::rowSds(value_samples))
  ) %>%
  select(-value, -value_samples) %>%
  arrange(output, longitude, latitude)

fst::write_fst(output, args$output)
