library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(ncdf4)
library(lubridate, warn.conflicts = FALSE)
library(parallel)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--oco2-observations')
parser$add_argument('--obspack-directory')
parser$add_argument('--tccon-sounding-directory')
parser$add_argument('--start-date')
parser$add_argument('--end-date')
parser$add_argument('--output')
args <- parser$parse_args()

log_info('Loading OCO-2 observations from {args$oco2_observations}')
with_nc_file(list(fn = args$oco2_observations), {
  v <- function(...) ncvar_get(fn, ...)
  oco2_soundings <- tibble(
    sounding_id = v('sounding_id'),
    assimilate = as.vector(v('assimilate_flag')),
    time = ncvar_get_time(fn, 'time'),
    longitude = as.vector(v('longitude')),
    latitude = as.vector(v('latitude')),
    value = as.vector(v('xco2')),
    measurement_error = as.vector(v('xco2_uncertainty')),
    model_error = as.vector(v('model_error')),
    oco2_dp = as.vector(v('dp')),
    oco2_log_dws = log(as.vector(v('dws'))),
    oco2_co2_grad_del = as.vector(v('co2_grad_del'))
  ) %>%
    mutate(
      oco2_operation_mode = c(
        'LN', 'LG', 'LT', 'LTT', 'ON', 'OG', 'OT', 'OTT'
      )[sounding_id %% 10]
    ) %>%
    mutate(
      oco2_operation_mode = factor(oco2_operation_mode)
    )
})

log_info('Loading ObsPack observations')
obspack_paths <- list.files(args$obspack_directory, full.names = TRUE)
obspack_times <- strptime(
  basename(obspack_paths),
  'obspack_co2_1_OCO2MIP_v3.2_2021-05-20.%Y%m%d.nc'
)
# NOTE(mgnb): add a few days buffer to account for overlap
obspack_paths <- obspack_paths[
  obspack_times >= (as.Date(args$start_date) - days(2))
  & obspack_times <= (as.Date(args$end_date) + days(2))
]
obspack_observations <- bind_rows(mclapply(obspack_paths, function(path) {
  log_trace('Opening {path}')
  fn <- nc_open(path)
  on.exit(nc_close(fn))
  v <- function(...) ncdf4::ncvar_get(fn, ...)
  models <- trimws(v('MIP_models'))
  tibble(
    obspack_id = as.vector(v('obspack_id')),
    assimilate = as.vector(v('CT_assim')),
    time = ncvar_get_time(fn, 'time'),
    longitude = as.vector(v('longitude')),
    latitude = as.vector(v('latitude')),
    altitude = as.vector(v('altitude')),
    value = 1e6 * as.vector(v('value')),
    # NOTE(mgnb): the MIP_MDM seems to have missing values
    model_error = 1e6 * as.vector(v('CT_MDM')),
    measurement_error = 0
  )
}, mc.cores = get_cores()))

log_info('Loading TCCON soundings')
tccon_paths <- list.files(args$tccon_sounding_directory, full.names = TRUE)
tccon_times <- strptime(basename(tccon_paths), 'tccon_timeavg_%Y%m%d.nc4')
tccon_paths <- tccon_paths[
  tccon_times >= args$start_date & tccon_times < args$end_date
]
tccon_soundings <- bind_rows(mclapply(tccon_paths, function(filename) {
  log_trace('Loading {filename}')
  with_nc_file(list(fn = filename), {
    v <- function(name, ...) ncvar_get(fn, sprintf('CO2/%s', name), collapse_degen = FALSE, ...)
    cdate <- v('cdate')
    tibble(
      observation_id = as.vector(v('obs_num')),
      time = ymd_hms(sprintf(
        '%04d-%02d-%02d %02d:%02d:%02d',
        cdate[1, ], cdate[2, ], cdate[3, ], cdate[4, ], cdate[5, ], cdate[6, ]
      )),
      longitude = as.vector(v('longitude')),
      latitude = as.vector(v('latitude')),
      value = as.vector(v('column_mixing'))
    )
  })
}, mc.cores = get_cores()))

log_info('Combining observations')
observations <- bind_rows(
  oco2_soundings %>%
    mutate(
      observation_id = as.character(sounding_id),
      observation_type = 'oco2'
    ) %>%
    select(-sounding_id),
  obspack_observations %>%
    mutate(
      observation_id = as.character(obspack_id),
      observation_type = 'obspack'
    ) %>%
    select(-obspack_id),
  tccon_soundings %>%
    mutate(
      observation_id = as.character(observation_id),
      observation_type = 'tccon',
      assimilate = 0
    )
) %>%
  select(observation_id, observation_type, time, everything()) %>%
  mutate(
    observation_type = factor(observation_type),
    overall_observation_mode = factor(case_when(
      observation_type == 'oco2'
        ~ as.character(oco2_operation_mode),
      observation_type == 'tccon'
        ~ 'TC',
      TRUE
        ~ 'IS'
    )),
    observation_group = factor(case_when(
      overall_observation_mode == 'IS'
        ~ stringr::str_split(observation_id, '~', simplify = TRUE)[, 2],
      overall_observation_mode %in% c('LN', 'LG')
        ~ '1_LNLG',
      overall_observation_mode == 'OG'
        ~ '2_OG',
      TRUE
        ~ as.character(overall_observation_mode)
    )),
    error = sqrt(
      model_error ^ 2 + measurement_error ^ 2
    ),
    type = stringr::str_split(
      as.character(observation_group),
      '_',
      simplify = TRUE
    )[, 3],
    hyperparameter_group = case_when(
      overall_observation_mode == 'IS'
        ~ stringr::str_split(type, '-', simplify = TRUE)[, 1],
      TRUE
        ~ as.character(observation_group)
    )
  ) %>%
  select(-type)

stopifnot(!any(is.na(observations$value)))
stopifnot(!any(is.na(observations$error[observations$assimilate == 1])))

log_info('Saving')
fst::write_fst(observations, args$output)

log_info('Done')
