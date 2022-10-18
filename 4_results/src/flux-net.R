library(argparse)
library(Matrix)
library(dplyr, warn.conflicts = FALSE)
library(Rcpp)
library(patchwork)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--area-1x1')
parser$add_argument('--perturbations-augmented')
parser$add_argument('--samples')
parser$add_argument('--region')
parser$add_argument('--output')
args <- parser$parse_args()

with_nc_file(list(fn = args$area_1x1), {
  longitude_area <- as.vector(ncdf4::ncvar_get(fn, 'lon'))
  latitude_area <- as.vector(ncdf4::ncvar_get(fn, 'lat'))
  area <- ncdf4::ncvar_get(fn, 'cell_area')
  area_1x1 <- expand.grid(
    longitude = longitude_area,
    latitude = latitude_area
  ) %>%
    mutate(area = as.vector(area))
})

if (!(args$region %in% names(REGION_PLOT_SETTINGS))) {
  stop('Invalid region')
}
plot_region <- REGION_PLOT_SETTINGS[[args$region]]
samples <- readRDS(args$samples)
perturbations_base <- fst::read_fst(args$perturbations_augmented)

area_495 <- (area_1x1 %>% filter(latitude == 49.5) %>% pull(area))[1]
area_505 <- (area_1x1 %>% filter(latitude == 50.5) %>% pull(area))[1]
area_both <- area_495 + area_505

# Splits grid cells that cross boundaries
perturbations_split <- bind_rows(
  perturbations_base %>%
    filter(latitude != 0),
  perturbations_base %>%
    filter(latitude == 0) %>%
    mutate(
      latitude = -0.5,
      latitude_bottom = -1,
      area = area / 2
    ),
  perturbations_base %>%
    filter(latitude == 0) %>%
    mutate(
      latitude = 0.5,
      latitude_bottom = 0,
      area = area / 2
    ),
  perturbations_base %>%
    filter(latitude == 50) %>%
    mutate(
      latitude = 49.5,
      latitude_bottom = 49,
      area = area * area_495 / area_both
    ),
  perturbations_base %>%
    filter(latitude == 50) %>%
    mutate(
      latitude = 50.5,
      latitude_bottom = 50,
      area = area * area_505 / area_both
    )
)

perturbations_region <- perturbations_split %>%
  filter(
    latitude_bottom >= plot_region$latitude_lower,
    latitude_bottom < plot_region$latitude_upper
  ) %>%
  group_by(inventory_minor_component_time, basis_vector) %>%
  summarise(value = KG_M2_S_TO_PGC_MONTH * sum(area * value)) %>%
  left_join(
    perturbations_split %>%
      distinct(inventory_minor_component_time, inventory, minor_component, time),
    by = 'inventory_minor_component_time'
  )

X_region <- with(perturbations_region, sparseMatrix(
  i = as.integer(inventory_minor_component_time),
  j = as.integer(basis_vector),
  x = value,
  dims = c(nlevels(inventory_minor_component_time), nlevels(basis_vector))
))

prior_emissions <- perturbations_region %>%
  group_by(inventory_minor_component_time, inventory, minor_component, time) %>%
  summarise(value = sum(value), .groups = 'drop') %>%
  select(-inventory_minor_component_time) %>%
  mutate(output = 'Bottom-up')

posterior_emissions <- prior_emissions %>%
  mutate(
    output = 'Posterior',
    value_prior = value,
    value = value_prior + as.vector(
      X_region[, as.integer(samples$alpha_df$basis_vector)]
      %*% samples$alpha_df$value
    ),
    value_samples = value_prior + as.matrix(
      X_region[, as.integer(samples$alpha_df$basis_vector)]
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
  mutate(
    inventory = factor(c(
      'bio_assim' = 'GPP',
      'bio_resp_tot' = 'Respiration',
      'nee' = 'NEE',
      'ocean' = 'Ocean'
    )[inventory], levels = c(
      'GPP',
      'Respiration',
      'NEE',
      'Ocean'
    )),
    minor_component = factor(c(
      'linear' = 'Linear',
      'periodic' = 'Seasonal',
      'residual' = 'Residual'
    )[as.character(minor_component)], levels = c(
      'Linear',
      'Seasonal',
      'Residual'
    ))
  )

emissions_net <- emissions %>%
  group_by(output, time, inventory) %>%
  summarise(
    value = sum(value),
    value_samples = t(colSums(value_samples)),
    .groups = 'drop'
  ) %>%
  mutate(
    value_q025 = matrixStats::rowQuantiles(value_samples, probs = 0.025),
    value_q975 = matrixStats::rowQuantiles(value_samples, probs = 0.975)
  )

output <- wrap_plots(lapply(sort(unique(emissions$inventory)), function(inventory_i) {
  ggplot(
    bind_rows(
      emissions_net %>%
        mutate(minor_component = 'Total'),
      emissions %>%
        filter(
          minor_component == 'Linear',
          !(inventory == 'Respiration' & output == 'Posterior'),
          !(inventory == 'Ocean' & output == 'Posterior')
        )
    ) %>%
      filter(inventory == inventory_i),
    aes(time)
  ) +
    geom_line(
      mapping = aes(
        y = value,
        colour = output,
        linetype = minor_component
      ),
      size = 0.4
    ) +
    geom_line(
      mapping = aes(
        y = value_q025,
        colour = output,
        linetype = minor_component
      ),
      size = 0.1
    ) +
    geom_line(
      mapping = aes(
        y = value_q975,
        colour = output,
        linetype = minor_component
      ),
      size = 0.1
    ) +
    scale_x_date(date_labels = '%Y-%m') +
    scale_colour_manual(values = c('black', '#ff4444')) +
    scale_fill_manual(values = c('black', '#ff4444')) +
    scale_linetype_manual(values = c('21', 'solid')) +
    labs(x = 'Time', y = 'Flux [PgC/month]', colour = NULL, fill = NULL, linetype = NULL) +
    guides(fill = 'none') +
    ggtitle(inventory_i)
}), ncol = 2, nrow = 2, guides = 'collect') &
  theme(
    plot.margin = margin(t = 1, r = 0, b = 0, l = 1, unit = 'mm'),
    plot.title = element_text(
      size = 11,
      margin = margin(0, 0, 5.5, 0, unit = 'points')
    ),
    legend.position = 'right',
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 10),
    strip.text = element_text(size = 8)
  )

output <- output +
  plot_annotation(
    title = sprintf('%s fluxes', plot_region$titlecase_title),
    theme = theme(
      plot.title = element_text(
        hjust = 0.4,
        size = 13
      )
    )
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$full_width,
  height = 8
)
