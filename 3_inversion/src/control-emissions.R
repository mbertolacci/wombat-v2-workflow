library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(ncdf4)
library(parallel)
library(lubridate, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--matched-runs')
parser$add_argument('--output')
args <- parser$parse_args()

log_debug('Loading emissions')
fn <- nc_open(file.path(args$matched_runs, 'base', 'monthly-fluxes.nc'))
on.exit(nc_close(fn))

v <- function(...) ncvar_get(fn, ...)

times <- as.Date(ncvar_get_time(fn, 'time'))
locations <- expand.grid(
  longitude = as.vector(v('lon')),
  latitude = as.vector(v('lat')),
  time = times
) %>%
  mutate(
    area = rep(
      as.vector(v('AREA')),
      length(times)
    ),
    cell_width = 2.5,
    cell_height = ifelse(
      abs(latitude) == 89.5,
      1,
      2
    )
  )

mapping <- data.frame(
  field = c(
    'EmisCO2_Biomass',
    'EmisCO2_Resp_Tot_Residual',
    'EmisCO2_Resp_Tot_Climatology',
    'EmisCO2_Assim_Residual',
    'EmisCO2_Assim_Climatology',
    'EmisCO2_Biofuel',
    'EmisCO2_Ocean_Residual',
    'EmisCO2_Ocean_Climatology',
    'EmisCO2_FossilFuel'
  ),
  inventory = factor(c(
    'biomass',
    'bio_resp_tot',
    'bio_resp_tot',
    'bio_assim',
    'bio_assim',
    'biofuel',
    'ocean',
    'ocean',
    'fossil_fuel'
  )),
  major_component = factor(c(
    'total',
    'residual',
    'climatology',
    'residual',
    'climatology',
    'total',
    'residual',
    'climatology',
    'total'
  ))
)

output <- bind_rows(lapply(mapping$field, function(field_i) {
  as_tibble(cbind(locations, data.frame(
    field = field_i,
    value = as.vector(v(field_i))
  )))
})) %>%
  left_join(mapping, by = 'field') %>%
  select(-field)

log_debug('Writing to {args$output}')
fst::write_fst(output, args$output)

log_debug('Done')
