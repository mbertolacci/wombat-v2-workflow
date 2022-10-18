library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(ncdf4)
library(tidyr)
library(parallel)

source(Sys.getenv('UTILS_PARTIAL'))

pluck <- function(x, name) lapply(x, getElement, name)

get_cell_width <- function(longitude) {
  rep(longitude[2] - longitude[1], length(longitude))
}

get_cell_height <- function(latitude) {
  unique_diff <- unique(diff(latitude))
  if (length(unique_diff) == 1) {
    rep(unique_diff, length(latitude))
  } else {
    # Half-height polar cells
    full_width <- max(unique_diff)
    c(full_width / 2, rep(full_width, length(latitude) - 2), full_width / 2)
  }
}

get_tccon_pressure_weights <- function(p) {
  # NOTE(mgnb): these are trapezoidal rule weights
  0.5 * abs(cbind(
    p[, 1] - p[, 2],
    matrixStats::rowDiffs(p, lag = 2),
    p[, ncol(p) - 1] - p[, ncol(p)]
  )) / p[, 1]
}

interpolate_vco2_pce <- function(from_vco2, from_pressure_edge, to_pressure_level) {
  n_levels_to <- ncol(to_pressure_level)

  k_for_l <- matrix(0L, nrow = nrow(from_vco2), ncol = n_levels_to)
  for (l in 1 : n_levels_to) {
    # NOTE(mgnb): if OCO-2 surface pressure is higher than GC, this is 0
    k_for_l[, l] <- rowSums(
      to_pressure_level[, l] <= from_pressure_edge
    )
  }
  to_vco2 <- to_pressure_level
  for (l in 1 : n_levels_to) {
    to_vco2[, l] <- from_vco2[cbind(
      seq_len(nrow(from_vco2)),
      # NOTE(mgnb): this extends the bottom GC layer down to the surface when
      # the to surface pressure is higher than the GC
      pmax(1, k_for_l[, l])
    )]
  }
  to_vco2
}

read_model <- function(species_conc_filename, meteorology_run) {
  species_conc_fn <- nc_open(species_conc_filename)

  variable_names <- names(species_conc_fn$var)
  species_names <- sort(sapply(
    strsplit(variable_names[startsWith(variable_names, 'SpeciesConc')], '_'),
    getElement,
    2
  ))
  time <- ncvar_get_time(species_conc_fn, 'time')
  time_width <- diff(time)[1]
  stopifnot(all(diff(time) == time_width))
  latitude <- as.vector(ncvar_get(species_conc_fn, 'lat'))
  cell_height <- get_cell_height(latitude)
  longitude <- as.vector(ncvar_get(species_conc_fn, 'lon'))
  cell_width <- get_cell_width(longitude)
  vertical_conc <- lapply(seq_along(species_names), function(j) {
    1e6 * ncvar_get(
      species_conc_fn,
      sprintf('SpeciesConc_%s', species_names[j])
    )
  })

  level_edge_filename <- file.path(
    meteorology_run,
    'output',
    stringr::str_replace(
      basename(species_conc_filename),
      'SpeciesConc',
      'LevelEdgeDiags'
    )
  )
  level_edge_fn <- nc_open(level_edge_filename)
  pressure_edge <- ncvar_get(level_edge_fn, 'Met_PEDGE')

  state_met_filename <- file.path(
    meteorology_run,
    'output',
    stringr::str_replace(
      basename(species_conc_filename),
      'SpeciesConc',
      'StateMet'
    )
  )
  state_met_fn <- nc_open(state_met_filename)
  box_height <- ncvar_get(state_met_fn, 'Met_BXHEIGHT')
  altitude <- aperm(
    apply(box_height, c(1, 2, 4), cumsum),
    c(2, 3, 1, 4)
  )

  nc_close(species_conc_fn)
  nc_close(level_edge_fn)
  nc_close(state_met_fn)

  list(
    time = time,
    time_width = time_width,
    latitude = latitude,
    cell_height = cell_height,
    longitude = longitude,
    cell_width = cell_width,
    species_names = species_names,
    vertical_conc = vertical_conc,
    pressure_edge = pressure_edge,
    altitude = altitude
  )
}

match_time <- function(time, model) {
  findInterval(
    time,
    c(model$time[1] - model$time_width / 2, model$time + model$time_width / 2),
    rightmost.closed = FALSE
  )
}

match_grid <- function(longitude, latitude, model) {
  latitude_index <- findInterval(
    latitude,
    c(model$latitude[1] - model$cell_height[1] / 2, model$latitude + model$cell_height / 2),
    rightmost.closed = TRUE
  )
  longitude_index <- findInterval(
    longitude,
    c(model$longitude[1] - model$cell_width[1] / 2, model$longitude + model$cell_width / 2),
    rightmost.closed = TRUE
  )
  longitude_index[longitude_index > length(model$longitude)] <- 1L

  data.frame(
    latitude = latitude_index,
    longitude = longitude_index
  )
}

match_soundings <- function(soundings, model) {
  period_start <- min(model$time) - model$time_width / 2
  period_end <- max(model$time) + model$time_width / 2

  log_trace(
    'Subsetting {soundings$observation_type[1]} soundings to between {period_start} and {period_end}'
  )
  soundings_part <- soundings %>%
    filter(
      time >= period_start,
      time < period_end
    )
  if (nrow(soundings_part) == 0) {
    return(tibble(
      observation_type = factor(),
      observation_id = factor(),
      model_time = POSIXct(),
      model_longitude = numeric(),
      model_latitude = numeric(),
      species = factor(),
      value = numeric()
    ))
  }

  log_trace('Computing match indices for {nrow(soundings_part)} soundings')
  indices <- match_grid(soundings_part$longitude, soundings_part$latitude, model)
  indices$time <- match_time(soundings_part$time, model)

  n_levels <- dim(model$vertical_conc[[1]])[3]

  log_trace('Subsetting vertical concentration')
  vco2_match <- array(NA, dim = c(
    length(indices$time),
    n_levels,
    length(model$species_names)
  ))
  for (j in seq_along(model$species_names)) {
    for (i in seq_along(indices$time)) {
      vco2_match[i, , j] <- model$vertical_conc[[j]][
        indices$longitude[i],
        indices$latitude[i],
        ,
        indices$time[i]
      ]
    }
  }
  log_trace('Subsetting pressure edges')
  pressure_edge_match <- matrix(
    NA,
    nrow = length(indices$time),
    ncol = n_levels + 1
  )
  for (i in seq_along(indices$time)) {
    pressure_edge_match[i, ] <- model$pressure_edge[
      indices$longitude[i],
      indices$latitude[i],
      ,
      indices$time[i]
    ]
  }

  output_base <- soundings_part %>%
    select(
      observation_type,
      observation_id = sounding_id
    ) %>%
    mutate(
      model_time = model$time[indices$time],
      model_longitude = model$longitude[indices$longitude],
      model_latitude = model$latitude[indices$latitude]
    )

  log_trace('Computing column concentration')
  bind_rows(lapply(seq_along(model$species_names), function(j) {
    vco2_match_j <- vco2_match[, , j]
    if (length(dim(vco2_match_j)) != 2) {
      vco2_match_j <- t(vco2_match_j)
    }
    vertical_conc_at_obs <- interpolate_vco2_pce(
      vco2_match_j,
      pressure_edge_match,
      soundings_part$pressure_levels
    )
    output_base$species <- model$species_names[j]
    output_base$value <- as.vector(soundings_part$xco2_apriori + rowSums(
      soundings_part$pressure_weights * soundings_part$averaging_kernel * (
        vertical_conc_at_obs - soundings_part$vco2_apriori
      ),
      na.rm = TRUE
    ))
    output_base
  })) %>%
    mutate(species = factor(species))
}

match_obspack <- function(observations, model) {
  empty_output <- tibble(
    observation_type = factor(),
    observation_id = factor(),
    observation_end_time = POSIXct(),
    model_time = POSIXct(),
    model_longitude = numeric(),
    model_latitude = numeric(),
    species = factor(),
    value = numeric(),
    sample_incomplete = logical(),
    n_samples = integer()
  )
  if (nrow(observations) == 0) {
    return(empty_output)
  }

  period_start <- min(model$time) - model$time_width / 2
  period_end <- max(model$time) + model$time_width / 2

  log_trace(
    'Subsetting {observations$observation_type[1]} observations to between {period_start} and {period_end}'
  )
  observations_part <- observations %>%
    filter(
      end_time >= period_start,
      start_time < period_end
    )
  if (nrow(observations_part) == 0) {
    return(empty_output)
  }

  log_trace('Computing match indices for {nrow(observations_part)} observations')
  n_times <- length(model$time)
  indices <- match_grid(observations_part$longitude, observations_part$latitude, model)
  central_time_index <- match_time(observations_part$time, model)
  start_time_index_raw <- match_time(observations_part$start_time, model)
  end_time_index_raw <- match_time(observations_part$end_time, model)
  start_time_index <- pmax(1L, start_time_index_raw)
  end_time_index <- pmin(n_times, end_time_index_raw)

  output_base <- observations_part %>%
    select(
      observation_type,
      observation_id = obspack_id,
      observation_altitude = altitude,
      observation_end_time = end_time
    ) %>%
    mutate(
      model_longitude = model$longitude[indices$longitude],
      model_latitude = model$latitude[indices$latitude],
      model_time = if_else(
        (central_time_index < 1) | (central_time_index > n_times),
        as_datetime(NA),
        model$time[pmax(1L, pmin(n_times, central_time_index))]
      ),
      sample_incomplete = (
        (start_time_index_raw < 1) | (end_time_index_raw > n_times)
      )
    )

  output_long <- output_base %>%
    select(observation_id, observation_altitude, observation_end_time) %>%
    mutate(
      longitude_index = indices$longitude,
      latitude_index = indices$latitude,
      time_indices = lapply(seq_along(start_time_index), function(i) {
        data.frame(time_index = start_time_index[i] : end_time_index[i])
      })
    ) %>%
    unnest(time_indices) %>%
    mutate(
      vertical_index = as.integer(sapply(seq_len(n()), function(i) {
        min(which(observation_altitude[i] < model$altitude[
          longitude_index[i],
          latitude_index[i],
          ,
          time_index[i]
        ]))
      }))
    ) %>%
    select(-observation_altitude)

  indices_4d <- output_long %>%
    select(
      longitude_index,
      latitude_index,
      vertical_index,
      time_index
    ) %>%
    as.matrix()

  long_factor <- factor(output_long$observation_id, output_base$observation_id)
  n_samples <- as.integer(table(long_factor))
  names(n_samples) <- NULL
  output_base$n_samples <- n_samples

  log_trace('Calculating mean concentrations')
  bind_rows(lapply(seq_along(model$species_names), function(j) {
    values <- as.vector(tapply(
      model$vertical_conc[[j]][indices_4d],
      long_factor,
      mean
    ))
    names(values) <- NULL
    output_base %>%
      mutate(
        species = model$species_names[j],
        value = values
      ) %>%
      select(-observation_altitude)
  })) %>%
    mutate(species = factor(species))
}

parser <- ArgumentParser()
parser$add_argument('--oco2-soundings')
parser$add_argument('--tccon-sounding-directory')
parser$add_argument('--obspack-directory')
parser$add_argument('--lauder')
parser$add_argument('--input-run')
parser$add_argument('--meteorology-run')
parser$add_argument('--output')
args <- parser$parse_args()

# args <- list(
#   oco2_soundings = 'data/OCO2_b10c_10sec_GOOD_r5.nc4',
#   tccon_sounding_directory = 'data/downloaded_20211217',
#   obspack_directory = 'data/obspack_co2_1_OCO2MIP_v3.2_2021-05-20/data/daily',
#   lauder = 'data/lauder_co2_2014_2021.50_ooofti_lhr.csv',
#   input_run = '1_transport/intermediates/runs-r10-r15-rNZ/residual_20210101_part001',
#   meteorology_run = '1_transport/intermediates/runs-r10-r15-rNZ/base',
#   output = '2_matching/intermediates/runs-r10-r15-rNZ/residual_20210101_part001'
# )

# args <- list(
#   oco2_soundings = 'data/OCO2_b10c_10sec_GOOD_r5.nc4',
#   tccon_sounding_directory = 'data/downloaded_20211217',
#   obspack_directory = 'data/obspack_co2_1_OCO2MIP_v3.2_2021-05-20/data/daily',
#   lauder = 'data/lauder_co2_2014_2021.50_ooofti_lhr.csv',
#   input_run = '1_transport/intermediates/runs-r10-r15-rNZ/residual_20181001_part001',
#   meteorology_run = '1_transport/intermediates/runs-r10-r15-rNZ/base',
#   output = '2_matching/intermediates/runs-r10-r15-rNZ/residual_20181001_part001'
# )

all_runs <- list.files(dirname(args$input_run))
input_runs <- sort(file.path(
  dirname(args$input_run),
  all_runs[startsWith(all_runs, basename(args$input_run))]
))

run_start_time <- NULL
run_end_time <- NULL
for (run in input_runs) {
  config_lines <- readLines(file.path(run, 'input.geos'), n = 10)
  run_start_time_current <- ymd_hms(config_lines[startsWith(config_lines, 'Start')])
  run_end_time_current <- ymd_hms(config_lines[startsWith(config_lines, 'End')])
  if (is.null(run_start_time) || run_start_time > run_start_time_current) {
    run_start_time <- run_start_time_current
  }
  if (is.null(run_end_time) || run_end_time < run_end_time_current) {
    run_end_time <- run_end_time_current
  }
}

log_info('Loading OCO-2 soundings from {args$oco2_soundings}')
with_nc_file(list(fn = args$oco2_soundings), {
  v <- function(...) ncvar_get(fn, ...)
  sigma_levels <- v('sigma_levels')

  oco2_soundings <- tibble(
    observation_type = factor('oco2'),
    sounding_id = factor(as.character(v('sounding_id'))),
    time = ncvar_get_time(fn, 'time'),
    longitude = as.vector(v('longitude')),
    latitude = as.vector(v('latitude')),
    pressure_surface = v('psurf'),
    pressure_weights = t(v('pressure_weight')),
    averaging_kernel = t(v('xco2_averaging_kernel')),
    xco2_apriori = v('xco2_apriori'),
    vco2_apriori = t(v('co2_profile_apriori')),
  ) %>%
    filter(time >= run_start_time, time < run_end_time)
})
oco2_soundings$pressure_levels <- (
  as.matrix(oco2_soundings$pressure_surface) %*% sigma_levels
)

log_info('Loading TCCON soundings')
tccon_paths <- list.files(args$tccon_sounding_directory, full.names = TRUE)
tccon_times <- strptime(basename(tccon_paths), 'tccon_timeavg_%Y%m%d.nc4')
tccon_paths <- tccon_paths[
  tccon_times >= run_start_time & tccon_times < run_end_time
]
tccon_soundings <- bind_rows(mclapply(tccon_paths, function(filename) {
  log_trace('Loading {filename}')
  with_nc_file(list(fn = filename), {
    v <- function(name, ...) ncvar_get(fn, sprintf('CO2/%s', name), collapse_degen = FALSE, ...)
    cdate <- v('cdate')
    tibble(
      observation_type = factor('tccon'),
      sounding_id = factor(as.character(as.vector(v('obs_num')))),
      time = ymd_hms(sprintf(
        '%04d-%02d-%02d %02d:%02d:%02d',
        cdate[1, ], cdate[2, ], cdate[3, ], cdate[4, ], cdate[5, ], cdate[6, ]
      )),
      longitude = as.vector(v('longitude')),
      latitude = as.vector(v('latitude')),
      pressure_levels = 1e-2 * t(v('p_levels_prior')),
      averaging_kernel = t(v('avg_kernel')),
      vh2o_wet_apriori = t(v('prior_h2o')),
      vco2_wet_apriori = t(v('prior_mixing'))
    )
  })
}, mc.cores = get_cores())) %>%
  mutate(pressure_surface = pressure_levels[, 1])
tccon_soundings$vco2_apriori <- with(tccon_soundings, {
  vco2_wet_apriori / (1 - vh2o_wet_apriori)
})
tccon_soundings$pressure_weights <- get_tccon_pressure_weights(tccon_soundings$pressure_levels)
tccon_soundings$xco2_apriori <- with(tccon_soundings, rowSums(
  pressure_weights * vco2_apriori
))

log_info('Loading ObsPack observations')
obspack_paths <- list.files(args$obspack_directory, full.names = TRUE)
obspack_times <- strptime(
  basename(obspack_paths),
  'obspack_co2_1_OCO2MIP_v3.2_2021-05-20.%Y%m%d.nc'
)
# NOTE(mgnb): add a few days buffer to account for overlap
obspack_paths <- obspack_paths[
  obspack_times >= (run_start_time - days(2))
  & obspack_times <= (run_end_time + days(2))
]
obspack_observations <- bind_rows(mclapply(obspack_paths, function(path) {
  log_trace('Opening {path}')
  fn <- nc_open(path)
  on.exit(nc_close(fn))
  v <- function(...) ncdf4::ncvar_get(fn, ...)
  tibble(
    observation_type = factor('obspack'),
    obspack_id = as.vector(v('obspack_id')),
    assimilate = as.vector(v('CT_assim')),
    time = ncvar_get_time(fn, 'time'),
    start_time = ncvar_get_time(fn, 'model_sample_window_start'),
    end_time = ncvar_get_time(fn, 'model_sample_window_end'),
    longitude = as.vector(v('longitude')),
    latitude = as.vector(v('latitude')),
    altitude = as.vector(v('altitude'))
  )
}, mc.cores = get_cores()))
if (nrow(obspack_observations) == 0) {
  obspack_observations <- tibble(
    observation_type = character(),
    obspack_id = character(),
    assimilate = numeric(),
    time = POSIXct(),
    start_time = POSIXct(),
    end_time = POSIXct(),
    longitude = numeric(),
    latitude = numeric(),
    altitude = numeric()
  )
}
obspack_observations <- obspack_observations %>%
  mutate(obspack_id = factor(obspack_id))
obspack_assimilate_split <- split(
  obspack_observations$obspack_id,
  obspack_observations$assimilate
)

log_info('Loading Lauder observations')
lauder_observations <- readr::read_csv(args$lauder, skip = 12, show_col_types = FALSE, lazy = FALSE) %>%
  mutate(
    observation_type = factor('lauder_is'),
    # Converts from NZST, which is +12 UTC
    time = ymd_hm(sprintf(
      '%04d-%02d-%02d %02d:30',
      YEAR,
      MON,
      DAY,
      `HOUR(START)`
    ), tz = 'UTC') - hours(12),
    assimilate = 3,
    start_time = time,
    end_time = time,
    altitude = 380,
    observation_id = factor(strftime(time, 'lauderis%Y%m%d%H%M', tz = 'UTC')),
    latitude = -45.038,
    longitude = 169.68
  ) %>%
  filter(time >= run_start_time, time < run_end_time) %>%
  select(
    observation_type,
    observation_id,
    assimilate,
    time,
    start_time,
    end_time,
    longitude,
    latitude,
    altitude
  )

for (resolution in c('hourly', 'daily')) {
  species_conc_files <- sort(do.call(c, lapply(input_runs, function(run) {
    list.files(
      file.path(run, 'output'),
      pattern = sprintf(
        'GEOSChem\\.SpeciesConc%s\\..*',
        paste0(toupper(substring(resolution, 1, 1)), substring(resolution, 2))
      ),
      full.names = TRUE
    )
  })))
  if (length(species_conc_files) == 0) {
    next
  }

  observation_dim <- ncdim_def('observation', '', 1L, unlim = TRUE, create_dimvar = FALSE)

  make_variable <- function(name, units, precision) {
    ncvar_def(name, units = units, dim = observation_dim, prec = precision, compression = 6)
  }

  base_output_vars <- list(
    make_variable('observation_type', 'factor', 'integer'),
    make_variable('observation_id', 'factor', 'integer'),
    make_variable('model_longitude', 'degrees_east', 'double'),
    make_variable('model_latitude', 'degrees_north', 'double'),
    make_variable('model_time', 'seconds since 1970-01-01 00:00:00 UTC', 'integer'),
    make_variable('species', 'factor', 'integer'),
    make_variable('value', 'ppm', 'double') #,
  )

  make_output_nc <- function(filename, df, id_name) {
    output <- nc_create(
      file.path(args$output, filename),
      base_output_vars,
      force_v4 = TRUE
    )
    attr(output, 'n_written') <- 0L
    output
  }

  oco2_fn <- make_output_nc(
    sprintf('oco2-%s.nc4', resolution),
    oco2_soundings,
    'sounding_id'
  )
  tccon_fn <- make_output_nc(
    sprintf('tccon-%s.nc4', resolution),
    tccon_soundings,
    'sounding_id'
  )
  lauder_fn <- make_output_nc(
    sprintf('lauder-%s.nc4', resolution),
    lauder_observations,
    'observation_id'
  )
  obspack_fns <- lapply(0 : 2, function(assimilate) {
    make_output_nc(
      sprintf('obspack-%s-assim-%d.nc4', resolution, assimilate),
      obspack_observations,
      'obspack_id'
    )
  })

  append_to_fn <- function(fn, df) {
    ncvar_put(
      fn,
      'observation_type',
      as.integer(df$observation_type),
      start = attr(fn, 'n_written') + 1L,
      count = nrow(df)
    )
    ncvar_put(
      fn,
      'observation_id',
      as.integer(df$observation_id),
      start = attr(fn, 'n_written') + 1L,
      count = nrow(df)
    )
    ncvar_put(
      fn,
      'model_longitude',
      df$model_longitude,
      start = attr(fn, 'n_written') + 1L,
      count = nrow(df)
    )
    ncvar_put(
      fn,
      'model_latitude',
      df$model_latitude,
      start = attr(fn, 'n_written') + 1L,
      count = nrow(df)
    )
    ncvar_put(
      fn,
      'model_time',
      as.integer(round(as.double(
        df$model_time - lubridate::ymd_hms('1970-01-01 00:00:00'),
        units = 'secs'
      ))),
      start = attr(fn, 'n_written') + 1L,
      count = nrow(df)
    )
    ncvar_put(
      fn,
      'species',
      as.integer(df$species),
      start = attr(fn, 'n_written') + 1L,
      count = nrow(df)
    )
    ncvar_put(
      fn,
      'value',
      df$value,
      start = attr(fn, 'n_written') + 1L,
      count = nrow(df)
    )
    attr(fn, 'n_written') <- attr(fn, 'n_written') + nrow(df)
    nc_sync(fn)
    fn
  }

  log_info('Processing {length(species_conc_files)} inputs for resolution = {resolution}')
  chunk_size <- 4 * get_cores()
  chunk_indices <- splitIndices(length(species_conc_files), length(species_conc_files) / chunk_size)
  obspack_buffer <- NULL
  for (chunk_index in seq_along(chunk_indices)) {
    chunk_indices_i <- chunk_indices[[chunk_index]]
    chunk_parts <- mclapply(species_conc_files[chunk_indices_i], function(species_conc_filename) {
      tryCatch({
        log_debug('[{round(pryr::mem_used() / 1024 ^ 2)} MB] Processing {species_conc_filename}')

        log_trace('Loading model variables')
        model <- read_model(species_conc_filename, args$meteorology_run)

        list(
          end_model_time = max(model$time) + model$time_width / 2,
          oco2_soundings = match_soundings(oco2_soundings, model),
          tccon_soundings = match_soundings(tccon_soundings, model),
          obspack = match_obspack(obspack_observations, model),
          lauder = match_obspack(
            lauder_observations %>% rename(obspack_id = observation_id),
            model
          )
        )
      }, error = function(e) {
        cat('Error in', species_conc_filename, '\n')
        print(e)
        stop(e)
      })
    }, mc.cores = get_cores())

    any_failed <- any(sapply(chunk_parts, is, 'try-error'))
    if (any_failed) {
      print(chunk_parts[sapply(chunk_parts, is, 'try-error')])
      stop('Error in matching')
    }

    log_trace('Appending to OCO-2 file')
    oco2_fn <<- append_to_fn(
      oco2_fn,
      bind_rows(pluck(chunk_parts, 'oco2_soundings'))
    )
    log_trace('Appending to TCCON file')
    tccon_fn <<- append_to_fn(
      tccon_fn,
      bind_rows(pluck(chunk_parts, 'tccon_soundings'))
    )
    log_trace('Appending to Lauder file')
    lauder_fn <<- append_to_fn(
      lauder_fn,
      bind_rows(pluck(chunk_parts, 'lauder'))
    )

    log_trace('Appending to ObsPack file')
    chunk_end_time <- chunk_parts[[length(chunk_parts)]]$end_model_time

    obspack_all <- c(obspack_buffer, pluck(chunk_parts, 'obspack'))

    obspack_complete <- lapply(obspack_all, function(part) {
      part %>%
        filter(!sample_incomplete)
    })

    obspack_ready_incomplete <- bind_rows(lapply(obspack_all, function(part) {
      part %>%
        filter(
          sample_incomplete,
          (chunk_index == length(chunk_indices))
          | (observation_end_time < chunk_end_time)
        )
    }))

    obspack_buffer <<- Filter(function(part) {
      nrow(part) > 0
    }, lapply(obspack_all, function(part) {
      part %>%
        filter(
          sample_incomplete,
          chunk_index != length(chunk_indices),
          observation_end_time >= chunk_end_time
        )
    }))

    obspack_ready_merged <- obspack_ready_incomplete %>%
      filter(!is.na(model_time)) %>%
      distinct(
        observation_id,
        species,
        .keep_all = TRUE
      ) %>%
      select(-n_samples, -value) %>%
      left_join(
        obspack_ready_incomplete %>%
          filter(!is.na(model_time)) %>%
          select(observation_id, species, n_samples, value) %>%
          group_by(observation_id, species) %>%
          summarise(
            value = sum(value * n_samples) / sum(n_samples),
            .groups = 'drop'
          ),
        by = c('observation_id', 'species')
      )

    obspack_to_sync <- c(obspack_complete, list(obspack_ready_merged))

    for (assimilate_plus_1 in 1L : 3L) {
      for (part in obspack_to_sync) {
        if (nrow(part) == 0) next
        obspack_fns[[assimilate_plus_1]] <- append_to_fn(
          obspack_fns[[assimilate_plus_1]],
          part[
            part$observation_id %in% obspack_assimilate_split[[assimilate_plus_1]],
          ]
        )
      }
    }

    print(gc())
  }

  nc_close(oco2_fn)
  nc_close(tccon_fn)
  nc_close(lauder_fn)
  for (assimilate_plus_1 in 1 : 3) {
    nc_close(obspack_fns[[assimilate_plus_1]])
  }

  species_names <- read_model(species_conc_files[1], args$meteorology_run)$species_names

  convert_nc_to_fst <- function(nc_filename, fst_filename, observation_df, id_name) {
    log_trace('Converting {nc_filename} to {fst_filename}')
    nc_path <- file.path(args$output, nc_filename)
    fn <- nc_open(nc_path)

    read_factor <- function(name, levels) {
      output <- as.integer(ncvar_get(fn, name))
      class(output) <- 'factor'
      attr(output, 'levels') <- levels
      output
    }
    read_double <- function(name) as.vector(ncvar_get(fn, name))

    fst::write_fst(tibble(
      observation_type = read_factor('observation_type', levels(observation_df$observation_type)),
      observation_id = read_factor('observation_id', levels(observation_df[[id_name]])),
      model_longitude = read_double('model_longitude'),
      model_latitude = read_double('model_latitude'),
      model_time = ncvar_get_time(fn, 'model_time'),
      species = read_factor('species', species_names),
      value = read_double('value')
    ), file.path(args$output, fst_filename))
    nc_close(fn)
    file.remove(nc_path)
  }

  convert_nc_to_fst(
    sprintf('oco2-%s.nc4', resolution),
    sprintf('oco2-%s.fst', resolution),
    oco2_soundings,
    'sounding_id'
  )

  convert_nc_to_fst(
    sprintf('tccon-%s.nc4', resolution),
    sprintf('tccon-%s.fst', resolution),
    tccon_soundings,
    'sounding_id'
  )

  convert_nc_to_fst(
    sprintf('lauder-%s.nc4', resolution),
    sprintf('lauder-%s.fst', resolution),
    lauder_observations,
    'observation_id'
  )

  for (assimilate in 0 : 2) {
    convert_nc_to_fst(
      sprintf('obspack-%s-assim-%s.nc4', resolution, assimilate),
      sprintf('obspack-%s-assim-%s.fst', resolution, assimilate),
      obspack_observations,
      'obspack_id'
    )
  }

  print(gc())
}

log_info('Done')
