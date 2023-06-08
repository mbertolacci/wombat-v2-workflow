library(argparse)
library(Matrix)
library(dplyr, warn.conflicts = FALSE)
library(patchwork)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--climatology-trend')
parser$add_argument('--six-year-average')
parser$add_argument('--region-sf')
parser$add_argument('--output')
args <- parser$parse_args()

climatology_trend <- fst::read_fst(args$climatology_trend)
six_year_average <- fst::read_fst(args$six_year_average)
region_sf <- readRDS(args$region_sf)

convert_to_sf <- function(df, mean_breaks, mean_limits, sd_breaks, sd_limits) {
  list(
    mean = df %>%
      filter(abs(latitude) != 89.5) %>%
      group_by(output) %>%
      group_map(~ {
        .x %>%
          mutate(
            value = discretise_by_breaks(value_mean, mean_breaks, mean_limits)
          ) %>%
          grid_df_to_sf('value') %>%
          mutate(
            output = .y$output
          )
      }) %>%
      bind_rows(),
    sd = df %>%
      filter(abs(latitude) != 89.5) %>%
      group_by(output) %>%
      group_map(~ {
        .x %>%
          mutate(
            value = discretise_by_breaks(value_sd, sd_breaks, sd_limits)
          ) %>%
          grid_df_to_sf('value') %>%
          mutate(
            output = .y$output
          )
      }) %>%
      bind_rows()
  )
}

average_mean_breaks <- round(seq(-0.6, 0.6, by = 0.1), 1)
average_mean_limits <- c(-0.6, 0.6)
average_sd_breaks <- round(seq(0, 0.1, by = 0.025), 3)
average_sd_limits <- c(0, 0.1 + 1.5 * 0.025)
six_year_average_sfs <- convert_to_sf(
  six_year_average,
  average_mean_breaks,
  average_mean_limits,
  average_sd_breaks,
  average_sd_limits
)

trend_mean_breaks <- seq(-0.02, 0.02, by = 0.005)
trend_mean_limits <- c(-0.02, 0.02) + c(-1, 1) * 0.005 * 1.5
trend_sd_breaks <- seq(0.001, 0.005, by = 0.001)
trend_sd_limits <- c(0, 0.005 + 1.5 * 0.001)
trend_sfs <- convert_to_sf(
  climatology_trend,
  trend_mean_breaks,
  trend_mean_limits,
  trend_sd_breaks,
  trend_sd_limits
)

average_bottom_up <- six_year_average_sfs$mean %>%
  filter(output == 'Bottom-up') %>%
  plot_map(value, average_mean_breaks, average_mean_limits, drop_second_labels = TRUE, show_excess = TRUE, label_precision = '1') +
  labs(
    fill = expression('Flux [kgCO'[2]*'/'*m^2*'/year]')
  ) +
  ggtitle('Bottom-up')

average_posterior_mean <- six_year_average_sfs$mean %>%
  filter(output == 'Posterior') %>%
  plot_map(value, average_mean_breaks, average_mean_limits, drop_second_labels = TRUE, show_excess = TRUE, label_precision = '1') +
  labs(
    fill = expression('Flux [kgCO'[2]*'/'*m^2*'/year]')
  ) +
  ggtitle('Posterior mean')

average_posterior_sd <- six_year_average_sfs$sd %>%
  filter(output == 'Posterior') %>%
  plot_map(value, average_sd_breaks, average_sd_limits, show_excess = TRUE, label_precision = '3', symmetric = FALSE) +
  labs(
    fill = expression('Posterior st. dev. of flux [kgCO'[2]*'/'*m^2*'/year]')
  ) +
  ggtitle('Posterior st. dev.')

trend_bottom_up <- trend_sfs$mean %>%
  filter(output == 'Bottom-up') %>%
  plot_map(value, trend_mean_breaks, trend_mean_limits, drop_second_labels = TRUE, show_excess = TRUE, label_precision = '3') +
  labs(
    fill = expression(beta['nee,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year'^2*']')
  ) +
  ggtitle('Bottom-up')

trend_posterior_mean <- trend_sfs$mean %>%
  filter(output == 'Posterior') %>%
  plot_map(value, trend_mean_breaks, trend_mean_limits, drop_second_labels = TRUE, show_excess = TRUE, label_precision = '3') +
  labs(
    fill = expression(beta['nee,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year'^2*']')
  ) +
  ggtitle('Posterior mean')

trend_posterior_sd <- trend_sfs$sd %>%
  filter(output == 'Posterior') %>%
  plot_map(value, trend_sd_breaks, trend_sd_limits, show_excess = TRUE, label_precision = '4', symmetric = FALSE) +
  labs(
    fill = expression('Posterior st. dev. of '*beta['nee,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year'^2*']')
  ) +
  ggtitle('Posterior st. dev.')

base_theme <- theme(
  legend.position = 'bottom',
  legend.margin = margin(t = -0.2, l = 0, b = -0.2, r = 0, unit = 'cm'),
  legend.title = element_text(size = 9),
  plot.title = element_text(
    hjust = 0.5,
    vjust = 1,
    size = 11 ,
    margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')
  ),
  plot.margin = margin(t = 0, b = 0.3, l = 0.1, r = 0, unit = 'cm')
)

top_theme <- base_theme +
  theme(
    legend.position = 'none',
    plot.margin = margin(t = 0.3, b = 0.2, l = 0.1, r = 0, unit = 'cm')
  )

sd_theme <- base_theme +
  theme(
    plot.margin = margin(t = 0.4, b = 0.2, l = 0.1, r = 0, unit = 'cm')
  )

average_column <- wrap_plots(
  average_bottom_up + top_theme,
  average_posterior_mean + base_theme,
  average_posterior_sd + sd_theme,
  ncol = 1
)

trend_column <- wrap_plots(
  trend_bottom_up + top_theme,
  trend_posterior_mean + base_theme,
  trend_posterior_sd + sd_theme,
  ncol = 1
)

output <- wrap_plots(
  wrap_plots(
    wrap_elements(
      panel = grid::textGrob('Average NEE flux from\nJanuary 2015 to December 2020')
    ),
    average_column,
    heights = c(0.08, 1)
  ),
  wrap_plots(
    wrap_elements(
      panel = grid::textGrob('Trend in NEE')
    ),
    trend_column,
    heights = c(0.08, 1)
  ),
  nrow = 1
)

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$full_width,
  height = 17.8
)
