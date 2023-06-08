source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))
library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(patchwork)
library(ncdf4)
library(tidyr)

parser <- ArgumentParser()
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

output <- example_cell %>%
  ungroup() %>%
  filter(
    latitude != 55.5
  ) %>%
  pivot_wider(
    names_from = 'major_component'
  ) %>%
  mutate(
    Total_Minus_Trend = Total - Linear
  ) %>%
  select(inventory, time, Total_Minus_Trend, Seasonal) %>%
  pivot_longer(-c(inventory, time)) %>%
  mutate(
    name = factor(c(
      'Total_Minus_Trend' = 'Total minus Linear',
      'Seasonal' = 'Seasonal'
    )[name], c(
      'Total minus Linear',
      'Seasonal'
    )),
    inventory = c(
      'bio_assim' = 'GPP',
      'bio_resp_tot' = 'Respiration'
    )[inventory]
  ) %>%
  ggplot(aes(time, value, colour = name)) +
    geom_line() +
    facet_wrap(~ inventory) +
    labs(x = 'Time', y = expression('Flux [kgC'*O[2]*'/'*m^2*'/year]'), colour = NULL) +
    theme(legend.position = 'bottom') +
    ggtitle('SiB4 Flux in 2.5° S 64.5° W')

log_info('Saving to {args$output}')
ggsave_fullwidth(
  args$output,
  output,
  height = 7
)

log_info('Done')
