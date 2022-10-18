library(argparse)
library(lubridate, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('CLIMATOLOGY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--year', nargs = '+', type = 'integer')
parser$add_argument('--split-year', action = 'store_true', default = FALSE)
parser$add_argument('--output', nargs = '+')
parser$add_argument('--climatology')
parser$add_argument('--field-name')
parser$add_argument('--time-resolution')
args <- parser$parse_args()

if (args$split_year) {
  for (year_index in seq_along(args$year)) {
    log_info('Expanding climatology for {args$year[year_index]}')
    if (args$time_resolution == 'hourly') {
      times <- seq(
        ymd_hm(sprintf('%d-01-01 00:00', args$year[year_index])),
        ymd_hm(sprintf('%d-12-31 23:00', args$year[year_index])),
        by = 'hour'
      )
    } else if (args$time_resolution == 'monthly') {
      times <- seq(
        ymd_hm(sprintf('%d-01-01 00:00', args$year[year_index])),
        ymd_hm(sprintf('%d-12-01 00:00', args$year[year_index])),
        by = 'month'
      )
    } else {
      stop('not supported')
    }

    expand_climatology(
      times,
      args$output[year_index],
      args$climatology,
      args$field_name,
      args$time_resolution,
      'start'
    )
  }
} else {
  if (args$time_resolution == 'hourly') {
    times <- seq(
      ymd_hm(sprintf('%d-01-01 00:00', head(args$year, 1))),
      ymd_hm(sprintf('%d-12-31 23:00', tail(args$year, 1))),
      by = 'hour'
    )
  } else if (args$time_resolution == 'monthly') {
    times <- seq(
      ymd_hm(sprintf('%d-01-01 00:00', head(args$year, 1))),
      ymd_hm(sprintf('%d-12-01 00:00', tail(args$year, 1))),
      by = 'month'
    )
  } else {
    stop('not supported')
  }

  expand_climatology(
    times,
    args$output,
    args$climatology,
    args$field_name,
    args$time_resolution,
    'start'
  )
}
