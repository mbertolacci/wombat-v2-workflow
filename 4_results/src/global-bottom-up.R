library(argparse)
library(dplyr)
library(patchwork)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--control-emissions')
parser$add_argument('--output')
args <- parser$parse_args()

control_emissions <- fst::read_fst(args$control_emissions)

bottom_up_emissions_global <- control_emissions %>%
  group_by(inventory, time) %>%
  summarise(
    value = KG_M2_S_TO_PGC_MONTH * sum(area * value),
    .groups = 'drop'
  ) %>%
  mutate(
    inventory = c(
      'biofuel' = 'biofuel',
      'biomass' = 'bioburn',
      'fossil_fuel' = 'fossil',
      'bio_assim' = 'gpp',
      'bio_resp_tot' = 'resp',
      'ocean' = 'ocean'
    )[as.character(inventory)]
  )

bottom_up_emissions_global <- bind_rows(
  bottom_up_emissions_global,
  bottom_up_emissions_global %>%
    filter(inventory %in% c('gpp', 'resp')) %>%
    group_by(time) %>%
    summarise(value = sum(value)) %>%
    mutate(inventory = 'nee (gpp + resp)')
)

bio_range <- bottom_up_emissions_global %>%
  filter(startsWith(as.character(inventory), 'bio')) %>%
  pull(value) %>%
  range()

output <- wrap_plots(lapply(c(
  'fossil', 'biofuel', 'bioburn', 'gpp', 'resp', 'nee (gpp + resp)', 'ocean'
), function(inventory_i) {
  output_i <- ggplot(
    bottom_up_emissions_global %>%
      filter(inventory == inventory_i),
    aes(time, value)
  ) +
    geom_line() +
    labs(x = 'Time', y = 'Flux [PgC/month]', colour = NULL, fill = NULL) +
    ggtitle(inventory_i) +
    theme(
      plot.margin = margin(t = 1, r = 3, b = 0, l = 3, unit = 'mm')
    )
  if (startsWith(inventory_i, 'bio')) {
    output_i <- output_i + ylim(bio_range)
  }
  output_i
}), ncol = 2, guides = 'collect')

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$full_width,
  height = 16
)
