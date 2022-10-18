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
parser$add_argument('--sib4-climatology-assim')
parser$add_argument('--sib4-climatology-resp-tot')
parser$add_argument('--sib4-monthly')
parser$add_argument('--output')
args <- parser$parse_args()

with_nc_file(list(fn = args$sib4_monthly), {
  assim <- ncvar_get(fn, 'assim')
  resp_tot <- ncvar_get(fn, 'resp_tot')
  latitude <- ncvar_get(fn, 'lat')
  longitude <- ncvar_get(fn, 'lon')
  time <- floor_date(
    ymd_hms('2000-01-01 00:00:00 UTC') + hours(24 * ncvar_get(fn, 'time')),
    'month'
  )
})

sib4_df <- data.frame(
  time = rep(time, each = length(latitude) * length(longitude)),
  latitude = rep(rep(latitude, each = length(longitude)), length(time)),
  longitude = rep(longitude, length(time) * length(latitude)),
  bio_assim = as.vector(assim),
  bio_resp_tot = as.vector(resp_tot)
) %>%
  pivot_longer(
    -c(time, longitude, latitude),
    names_to = 'inventory'
  )

example_cell <- data.frame(
  longitude = c(95.5, -64.5),
  latitude = c(55.5, -2.5)
) %>%
  left_join(
    sib4_df,
    by = c('longitude', 'latitude')
  ) %>%
  group_by(longitude, latitude, inventory) %>%
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
  args$sib4_climatology_assim,
  'coefficient',
  include_time = FALSE,
  include_variable = TRUE
), 'origin'))

sib4_climatology <- bind_rows(
  read_climatology(args$sib4_climatology_assim) %>%
    mutate(inventory = 'bio_assim'),
  read_climatology(args$sib4_climatology_resp_tot) %>%
    mutate(inventory = 'bio_resp_tot')
)

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

sib4_global_total <- sib4_df %>%
  left_join(area_df, by = c('longitude', 'latitude')) %>%
  group_by(time, inventory) %>%
  summarise(
    value = KG_M2_S_TO_PGC_MONTH * sum(value * area)
  )

sib4_global_climatology <- sib4_climatology %>%
  left_join(area_df, by = c('longitude', 'latitude')) %>%
  group_by(variable, inventory) %>%
  summarise(
    value = KG_M2_S_TO_PGC_MONTH * sum(value * area)
  )

sib4_global <- sib4_global_climatology %>%
  group_by(inventory) %>%
  group_modify(~ {
    times <- sort(unique(sib4_global_total$time))
    X <- climatology_design_matrix(
      time = times,
      origin = origin,
      time_alignment = 'start',
      time_resolution = 'monthly',
      harmonics = 3
    )
    data.frame(
      time = times,
      Linear = as.vector(X[, 1 : 2] %*% .x$value[1 : 2]),
      Seasonal = as.vector(X[, 3 : ncol(X)] %*% .x$value[3 : ncol(X)])
    )
  }) %>%
  left_join(
    sib4_global_total %>%
      rename(Total = value),
    by = c('inventory', 'time')
  ) %>%
  mutate(
    Residual = Total - Linear - Seasonal
  ) %>%
  pivot_longer(
    -c(inventory, time),
    names_to = 'major_component'
  )

plot_series <- function(df) {
  input_df <- df %>%
    mutate(
      inventory = c(
        'bio_assim' = 'GPP',
        'bio_resp_tot' = 'Respiration'
      )[inventory],
      major_component = factor(
        major_component,
        levels = c('Total', 'Linear', 'Seasonal', 'Residual')
      )
    )

  modify_major_component <- function(x) {
    factor(x, levels = levels(input_df$major_component))
  }

  ggplot(input_df, aes(time, value, colour = inventory)) +
    geom_rect(
      data = data.frame(
        start = ymd_hms('2015-01-01 00:00:00'),
        end = ymd_hms('2020-12-31 00:00:00')
      ),
      mapping = aes(
        xmin = start,
        xmax = end,
        ymin = -Inf,
        ymax = Inf
      ),
      colour = NA,
      fill = 'black',
      alpha = 0.1,
      inherit.aes = FALSE
    ) +
    geom_line(size = 0.4) +
    facet_grid(major_component ~ ., scales = 'free_y') +
    labs(x = 'Time', colour = 'Component', y = NULL) +
    scale_colour_manual(
      values = c('#35978f', '#bf812d')
    )
}

output1 <- example_cell %>%
  filter(latitude == 55.5) %>%
  plot_series() +
    labs(y = expression('Flux [kgC'*O[2]*'/'*m^2*'/year]')) +
    ggtitle('55.5° N 3.5° W')

output2 <- example_cell %>%
  filter(latitude != 55.5) %>%
  plot_series() +
    labs(y = NULL) +
    ggtitle('2.5° S 64.5° W')

output3 <- sib4_global %>%
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
    axis.text.x = element_text(size = 10, angle = 33, hjust = 1),
    axis.text.y = element_text(size = 9),
    legend.position = 'bottom',
    legend.margin = margin(t = -3, r = 0, b = 0, l = 0, unit = 'mm')
  )

log_info('Saving to {args$output}')
ggsave_fullwidth(
  args$output,
  output,
  height = 8.8
)

log_info('Done')
