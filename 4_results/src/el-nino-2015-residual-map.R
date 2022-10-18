library(argparse)
library(Matrix)
library(dplyr, warn.conflicts = FALSE)
library(patchwork)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations-augmented')
parser$add_argument('--region-sf')
parser$add_argument('--samples')
parser$add_argument('--output')
args <- parser$parse_args()

samples <- readRDS(args$samples)
region_sf <- readRDS(args$region_sf)
control_emissions <- fst::read_fst(args$control_emissions)
perturbations <- fst::read_fst(args$perturbations_augmented)

perturbations_base <- perturbations %>%
  filter(
    inventory != 'ocean',
    component == 'residual',
    time >= '2015-12-01',
    time < '2016-03-01'
  ) %>%
  mutate(value = PER_SECONDS_TO_PER_YEAR * value)

el_nino_residual_parts <- bind_rows(
  perturbations_base %>%
    mutate(output = 'Bottom-up'),
  perturbations_base %>%
    left_join(
      samples$alpha_df %>%
        select(basis_vector, alpha = value),
      by = 'basis_vector'
    ) %>%
    mutate(
      output = 'Posterior mean',
      value = (1 + alpha) * value
    )
) %>%
  group_by(output, inventory, longitude, latitude, time) %>%
  summarise(
    value = sum(value)
  ) %>%
  group_by(output, inventory, longitude, latitude) %>%
  summarise(
    value = mean(value),
    .groups = 'drop'
  )

el_nino_residual <- bind_rows(
  el_nino_residual_parts,
  el_nino_residual_parts %>%
    group_by(output, longitude, latitude) %>%
    summarise(
      value = sum(value)
    ) %>%
    mutate(inventory = 'nee')
) %>%
  mutate(
    inventory = factor(c(
      'bio_assim' = 'GPP',
      'bio_resp_tot' = 'Respiration',
      'nee' = 'NEE'
    )[inventory], c(
      'GPP',
      'Respiration',
      'NEE'
    ))
  )

flux_breaks <- seq(-4.5, 4.5, by = 0.5)
flux_limits <- c(-4.5, 4.5)

el_nino_residual_mean_sf <- el_nino_residual %>%
  group_by(inventory, output) %>%
  group_map(~ {
    control_emissions %>%
      distinct(longitude, latitude) %>%
      arrange(longitude, latitude) %>%
      filter(abs(latitude) != 89.5) %>%
      left_join(.x, by = c('longitude', 'latitude')) %>%
      mutate(
        value = discretise_by_breaks(value, flux_breaks, flux_limits)
      ) %>%
      grid_df_to_sf('value') %>%
      mutate(
        inventory = .y$inventory,
        output = .y$output
      )
  }) %>%
  bind_rows()

output <- plot_map(
  el_nino_residual_mean_sf,
  value,
  flux_breaks,
  flux_limits,
  show_excess = FALSE,
  label_precision = 1,
  bar_width = 16,
  drop_second_labels = TRUE
) +
  facet_grid(inventory ~ output) +
  labs(
    fill = expression('Flux [kgCO'[2]*'/'*m^2*'/year]')
  ) +
  ggtitle('Average residual flux over December 2015 to February 2016') +
  theme(
    legend.position = 'bottom',
    plot.title = element_text(size = 13, hjust = 0.5),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(size = 12)
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 16
)
