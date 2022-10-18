library(argparse)
library(Matrix)
library(dplyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations-augmented')
parser$add_argument('--samples')
parser$add_argument('--output')
args <- parser$parse_args()

samples <- readRDS(args$samples)
control_emissions <- fst::read_fst(args$control_emissions)
perturbations <- fst::read_fst(args$perturbations_augmented)

perturbations_base <- perturbations %>%
  filter(
    time >= '2015-01-01',
    time < '2021-01-01',
    inventory %in% c('bio_assim', 'bio_resp_tot')
  ) %>%
  mutate(
    value = PER_SECONDS_TO_PER_YEAR * value
  )

locations <- control_emissions %>%
  distinct(longitude, latitude) %>%
  arrange(longitude, latitude) %>%
  mutate(location_index = 1 : n())

X_alpha_to_six_year_mean <- with(
  perturbations_base %>%
    group_by(basis_vector, longitude, latitude) %>%
    summarise(value = sum(value) / 72, .groups = 'drop') %>%
    left_join(locations, by = c('longitude', 'latitude')),
  sparseMatrix(
    i = location_index,
    j = as.integer(basis_vector),
    x = value,
    dims = c(nrow(locations), nlevels(basis_vector))
  )
)[, as.integer(samples$alpha_df$basis_vector)]

six_year_average_bottom_up <- locations %>%
  select(-location_index) %>%
  left_join(
    perturbations_base %>%
      group_by(longitude, latitude, time) %>%
      summarise(
        value = sum(value),
        .groups = 'drop'
      ) %>%
      group_by(longitude, latitude) %>%
      summarise(
        value = mean(value),
        .groups = 'drop'
      ),
    by = c('longitude', 'latitude')
  )

alpha_mean <- samples$alpha_df$value
alpha_samples <- samples$alpha_df$value_samples[
  ,
  floor(seq(1, ncol(samples$alpha_df$value_samples), length.out = 100))
]
value_tilde_mean <- as.vector(X_alpha_to_six_year_mean %*% alpha_mean)
value_tilde_samples <- as.matrix(X_alpha_to_six_year_mean %*% alpha_samples)

output <- bind_rows(
  six_year_average_bottom_up %>%
    mutate(output = 'Bottom-up'),
  six_year_average_bottom_up %>%
    mutate(
      output = 'Posterior',
      value_samples = value + value_tilde_samples,
      value = value + value_tilde_mean,
      value_sd = matrixStats::rowSds(value_samples)
    )
) %>%
  select(-value_samples) %>%
  rename(value_mean = value)

fst::write_fst(output, args$output)
