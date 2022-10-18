library(argparse)
library(Matrix)
library(dplyr, warn.conflicts = FALSE)
library(Rcpp)
library(patchwork)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--basis-vectors')
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations')
parser$add_argument('--samples')
parser$add_argument('--output')
args <- parser$parse_args()

basis_vectors <- fst::read_fst(args$basis_vectors)
samples <- readRDS(args$samples)
control_emissions <- fst::read_fst(args$control_emissions)
perturbations <- fst::read_fst(args$perturbations)
cell_area <- control_emissions %>%
  distinct(longitude, latitude, area)

perturbations_base <- perturbations %>%
  add_basis_vector(basis_vectors) %>%
  mutate(
    minor_component = factor(
      case_when(
        component == 'residual' ~ 'residual',
        component %in% c('intercept', 'trend') ~ 'linear',
        TRUE ~ 'periodic'
      ),
      c('linear', 'periodic', 'residual')
    ),
    inventory_minor_component_time = interaction(
      inventory,
      minor_component,
      time,
      drop = TRUE
    )
  ) %>%
  left_join(cell_area, by = c('longitude', 'latitude'))

perturbations_global <- perturbations_base %>%
  group_by(inventory_minor_component_time, basis_vector) %>%
  summarise(value = KG_M2_S_TO_PGC_MONTH * sum(area * value)) %>%
  left_join(
    perturbations_base %>%
      distinct(inventory_minor_component_time, inventory, minor_component, time),
    by = 'inventory_minor_component_time'
  )

X_global <- with(perturbations_global, sparseMatrix(
  i = as.integer(inventory_minor_component_time),
  j = as.integer(basis_vector),
  x = value,
  dims = c(nlevels(inventory_minor_component_time), nlevels(basis_vector))
))

prior_emissions <- perturbations_global %>%
  group_by(inventory_minor_component_time, inventory, minor_component, time) %>%
  summarise(value = sum(value), .groups = 'drop') %>%
  select(-inventory_minor_component_time) %>%
  mutate(output = 'Bottom-up\n(prior mode)')

posterior_emissions <- prior_emissions %>%
  mutate(
    output = 'Posterior',
    value_prior = value,
    value = value_prior + as.vector(
      X_global[, as.integer(samples$alpha_df$basis_vector)]
      %*% samples$alpha_df$value
    ),
    value_samples = value_prior + as.matrix(
      X_global[, as.integer(samples$alpha_df$basis_vector)]
      %*% samples$alpha_df$value_samples
    ),
    value_q025 = matrixStats::rowQuantiles(value_samples, probs = 0.025),
    value_q975 = matrixStats::rowQuantiles(value_samples, probs = 0.975)
  ) %>%
  select(-value_prior)

emissions <- bind_rows(
  prior_emissions,
  posterior_emissions
) %>%
  {
    x <- .

    bind_rows(
      x,
      x %>%
        filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
        group_by(output, time, minor_component) %>%
        summarise(
          value = sum(value),
          value_samples = t(colSums(value_samples)),
          .groups = 'drop'
        ) %>%
        mutate(
          inventory = 'nee',
          value_q025 = matrixStats::rowQuantiles(value_samples, probs = 0.025),
          value_q975 = matrixStats::rowQuantiles(value_samples, probs = 0.975)
        )
    )
  } %>%
  filter(!(
    (inventory == 'bio_resp_tot' & minor_component == 'linear' & output == 'Posterior')
    | (inventory == 'ocean' & minor_component %in% c('linear', 'periodic') & output == 'Posterior')
  )) %>%
  mutate(
    inventory = factor(c(
      'bio_assim' = 'bioassim (GPP)',
      'bio_resp_tot' = 'bioresp',
      'nee' = 'nee',
      'ocean' = 'ocean'
    )[inventory], levels = c(
      'bioassim (GPP)',
      'bioresp',
      'nee',
      'ocean'
    ))
  )

output <- wrap_plots(lapply(sort(unique(emissions$inventory)), function(inventory_i) {
  ggplot(
    emissions %>%
      filter(inventory == inventory_i),
    aes(time)
  ) +
    geom_line(
      mapping = aes(
        y = value,
        colour = output,
        linetype = output
      ),
      size = 0.4
    ) +
    geom_ribbon(
      mapping = aes(
        ymin = value_q025,
        ymax = value_q975,
        fill = output
      ),
      alpha = 0.3
    ) +
    facet_grid(minor_component ~ ., scales = 'free_y') +
    scale_colour_manual(values = c('black', '#ff4444')) +
    scale_fill_manual(values = c('black', '#ff4444')) +
    scale_linetype_manual(values = c('41', 'solid')) +
    labs(x = 'Time', y = 'Flux [PgC/month]', colour = NULL, fill = NULL, linetype = NULL) +
    guides(fill = 'none') +
    ggtitle(inventory_i)
}), ncol = 2, nrow = 2, guides = 'collect') &
  theme(
    plot.margin = margin(t = 1, r = 0, b = 0, l = 1, unit = 'mm'),
    legend.position = 'right',
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 10),
    strip.text = element_text(size = 9)
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$full_width,
  height = 13
)
