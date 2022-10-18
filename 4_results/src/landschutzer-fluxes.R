source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))
library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(patchwork)
library(ncdf4)
library(tidyr)

parser <- ArgumentParser()
parser$add_argument('--area-1x1')
parser$add_argument('--landschutzer')
parser$add_argument('--landschutzer-climatology')
parser$add_argument('--show-full', action = 'store_true', default = FALSE)
parser$add_argument('--output')
args <- parser$parse_args()

with_nc_file(list(fn = args$landschutzer), {
  fgco2_raw <- ncvar_get(fn, 'fgco2_raw')
  latitude <- ncvar_get(fn, 'lat')
  longitude <- ncvar_get(fn, 'lon')
  time <- ymd_hms('2000-01-01 00:00:00 UTC') + seconds(ncvar_get(fn, 'time'))
})

landschutzer_df <- data.frame(
  time = rep(time, each = length(latitude) * length(longitude)),
  latitude = rep(rep(latitude, each = length(longitude)), length(time)),
  longitude = rep(longitude, length(time) * length(latitude)),
  value = as.vector(fgco2_raw)
)

example_cell <- data.frame(
  longitude = c(-90.5, -28.5),
  latitude = c(-10.5, 37.5)
) %>%
  left_join(
    landschutzer_df,
    by = c('longitude', 'latitude')
  ) %>%
  group_by(longitude, latitude) %>%
  group_modify(~ {
    X <- climatology_design_matrix(
      time = .x$time,
      origin = .x$time[1],
      time_alignment = 'start',
      time_resolution = 'monthly',
      harmonics = 3
    )
    beta_hat <- as.vector(solve(crossprod(X), crossprod(X, .x$value)))
    bind_rows(
      data.frame(
        time = .x$time,
        major_component = 'Total',
        value = .x$value
      ),
      data.frame(
        time = .x$time,
        major_component = 'Linear',
        value = as.vector(X[, 1 : 2] %*% beta_hat[1 : 2])
      ),
      data.frame(
        time = .x$time,
        major_component = 'Seasonal',
        value = as.vector(X[, 3 : ncol(X)] %*% beta_hat[3 : ncol(X)])
      ),
      data.frame(
        time = .x$time,
        major_component = 'Residual',
        value = .x$value - as.vector(X %*% beta_hat)
      )
    )
  }) %>%
  mutate(
    # Convert kg/m^2/s to kg/m^2/year
    value = 365.24 * 24 * 60 * 60 * value
  )

origin <- ymd_hms(attr(read_gridded_data(
  args$landschutzer_climatology,
  'coefficient',
  include_time = FALSE,
  include_variable = TRUE
), 'origin'))

landschutzer_climatology <- read_climatology(args$landschutzer_climatology)

with_nc_file(list(fn = args$area_1x1), {
  longitude_area <- as.vector(ncdf4::ncvar_get(fn, 'lon'))
  latitude_area <- as.vector(ncdf4::ncvar_get(fn, 'lat'))
  area <- ncdf4::ncvar_get(fn, 'cell_area')
})

area_df <- expand.grid(
  longitude = longitude_area,
  latitude = latitude_area
) %>%
  mutate(area = as.vector(area))

KG_M2_S_TO_PGC_MONTH <- 12.01 / 44.01 * 1e-12 * 365.25 * 24 * 60 * 60 / 12

landschutzer_global_total <- landschutzer_df %>%
  left_join(area_df, by = c('longitude', 'latitude')) %>%
  group_by(time) %>%
  summarise(
    value = KG_M2_S_TO_PGC_MONTH * sum(value * area)
  )

landschutzer_global_climatology <- landschutzer_climatology %>%
  left_join(area_df, by = c('longitude', 'latitude')) %>%
  group_by(variable) %>%
  summarise(
    value = KG_M2_S_TO_PGC_MONTH * sum(value * area)
  )

times <- sort(unique(landschutzer_global_total$time))
X <- climatology_design_matrix(
  time = times,
  origin = origin,
  time_alignment = 'start',
  time_resolution = 'monthly',
  harmonics = 2
)
landschutzer_total <- data.frame(
  time = times,
  Linear = as.vector(X[, 1 : 2] %*% landschutzer_global_climatology$value[1 : 2]),
  Seasonal = as.vector(X[, 3 : ncol(X)] %*% landschutzer_global_climatology$value[3 : ncol(X)])
) %>%
  left_join(
    landschutzer_global_total %>%
      rename(Total = value),
    by = 'time'
  ) %>%
  mutate(
    Residual = Total - Linear - Seasonal
  ) %>%
  pivot_longer(
    -time,
    names_to = 'major_component'
  )

plot_series <- function(df) {
  input_df <- df %>%
    mutate(
      major_component = factor(
        major_component,
        levels = c('Total', 'Linear', 'Seasonal', 'Residual')
      )
    )
  if (!args$show_full) {
    input_df <- input_df %>% filter(time >= '2010-01-01')
  }

  modify_major_component <- function(x) {
    factor(x, levels = levels(input_df$major_component))
  }

  ggplot(input_df, aes(time, value)) +
    geom_line(size = 0.4) +
    facet_grid(major_component ~ ., scales = 'free_y') +
    labs(x = 'Time', y = NULL)
}

output1 <- example_cell %>%
  filter(latitude == -10.5) %>%
  plot_series() +
    labs(y = expression('Flux [kgC'*O[2]*'/'*m^2*'/year]')) +
    ggtitle('10.5° S 90.5° W')

output2 <- example_cell %>%
  filter(latitude != -10.5) %>%
  plot_series() +
    labs(y = NULL) +
    ggtitle('37.5° N 28.5° W')

output3 <- landschutzer_total %>%
  plot_series() +
    labs(y = 'Flux [PgC/month]') +
    ggtitle('Global total')

output <- wrap_plots(
  output1,
  output2,
  output3,
  nrow = 1,
  guides = 'collect'
) &
  theme(
    legend.position = 'bottom',
    legend.margin = margin(t = 0, unit = 'mm')
  )

log_info('Saving to {args$output}')
ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 8.5
)

log_info('Done')
