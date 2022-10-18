estimate_climatology <- function(
  input_filename,
  field_name,
  time_resolution = c('hourly', 'monthly'),
  time_alignment = c('start', 'middle', 'end'),
  origin,
  harmonics,
  output_filename
) {
  time_resolution <- match.arg(time_resolution)
  time_alignment <- match.arg(time_alignment)

  log_debug('Loading {field_name} from {input_filename}')
  input <- read_gridded_data(input_filename, field_name)

  log_debug('Decomposing climatology')
  outputs <- array(
    0,
    dim = c(length(input$longitude), length(input$latitude), 4 * harmonics + 2)
  )
  dimnames(outputs) <- list(
    longitude = NULL,
    latitude = NULL,
    variable = c(
      'intercept',
      'trend',
      as.vector(sapply(seq_len(harmonics), function(k) {
        sprintf(c('cos12_%d_intercept', 'sin12_%d_intercept'), k)
      })),
      as.vector(sapply(seq_len(harmonics), function(k) {
        sprintf(c('cos12_%d_trend', 'sin12_%d_trend'), k)
      }))
    )
  )

  X <- climatology_design_matrix(
    time = input$time,
    origin = origin,
    time_alignment = time_alignment,
    time_resolution = time_resolution,
    harmonics = harmonics
  )
  for (latitude_index in seq_along(input$latitude)) {
    for (longitude_index in seq_along(input$longitude)) {
      y <- input$value[longitude_index, latitude_index, ]
      beta_hat <- as.vector(solve(crossprod(X), crossprod(X, y)))
      outputs[longitude_index, latitude_index, ] <- beta_hat
    }
  }

  log_debug('Writing to {output_filename}')
  dimensions <- lapply(list(
    list(
      name = 'lon',
      units = 'degrees_east',
      longname = 'Longitude',
      vals = input$longitude
    ),
    list(
      name = 'lat',
      units = 'degrees_north',
      longname = 'Latitude',
      vals = input$latitude
    ),
    list(
      name = 'variable',
      units = '',
      longname = 'Variable',
      vals = seq_len(length(dimnames(outputs)$variable)),
      create_dimvar = FALSE
    )
  ), function(args) do.call(ncdf4::ncdim_def, args))

  coefficient_variable <- ncdf4::ncvar_def(
    'coefficient',
    units = input$value_units,
    longname = 'Coefficients',
    dim = dimensions,
    prec = 'double',
    compression = 6
  )
  variable_variable <- ncdf4::ncvar_def(
    'variable',
    units = '',
    longname = 'Variables',
    dim = list(
      ncdf4::ncdim_def(
        'nchar',
        '',
        seq_len(max(nchar(dimnames(outputs)$variable))),
        create_dimvar = FALSE
      ),
      dimensions[[3]]
    ),
    prec = 'char',
    compression = 6
  )

  output_fn <- ncdf4::nc_create(
    output_filename,
    list(coefficient_variable, variable_variable),
    force_v4 = TRUE
  )
  on.exit(ncdf4::nc_close(output_fn))
  ncdf4::ncvar_put(output_fn, coefficient_variable, outputs)
  ncdf4::ncvar_put(output_fn, variable_variable, dimnames(outputs)$variable)
  ncdf4::ncatt_put(output_fn, 0, 'origin', format(origin, '%Y-%m-%d %H:%M:%S'))
  ncdf4::ncatt_put(output_fn, 0, 'harmonics', harmonics)

  invisible(NULL)
}

decompose_climatology <- function(
  input_filenames,
  output_filenames,
  climatology_filename,
  field_name,
  time_resolution = c('hourly', 'monthly'),
  time_alignment = c('start', 'middle', 'end'),
  fields
) {
  has_fields <- !missing(fields)
  time_resolution <- match.arg(time_resolution)
  time_alignment <- match.arg(time_alignment)

  log_debug('Loading climatology')
  climatology <- read_gridded_data(
    climatology_filename,
    'coefficient',
    include_time = FALSE,
    include_variable = TRUE
  )

  for (file_index in seq_along(input_filenames)) {
    log_debug('Processing {input_filenames[file_index]}')
    input <- read_gridded_data(input_filenames[file_index], field_name)
    with_nc_file(list(fn = input_filenames[file_index]), {
      dimensions <- fn$var[[field_name]]$dim
    })

    log_debug('Calculating design matrix over {length(input$time)} time periods')
    X <- climatology_design_matrix(
      time = input$time,
      origin = ymd_hms(attr(climatology, 'origin')),
      time_alignment = time_alignment,
      time_resolution = time_resolution,
      harmonics = attr(climatology, 'harmonics')
    )

    lapply_variable <- function(fn) {
      lapply(seq_along(climatology$variable), fn)
    }

    log_trace('Defining variables')
    variables <- lapply_variable(function(variable_index) {
      variable_name <- climatology$variable[variable_index]
      ncdf4::ncvar_def(
        variable_name,
        units = input$value_units,
        longname = sprintf('%s %s', input$value_longname, variable_name),
        dim = dimensions,
        prec = 'double',
        compression = 6,
        chunksizes = c(dim(climatology$value)[1 : 2], 4)
      )
    })
    variables <- Filter(Negate(is.null), variables)
    residual_variable <- ncdf4::ncvar_def(
      'residual',
      units = input$value_units,
      longname = sprintf('%s residual', input$value_longname),
      dim = dimensions,
      prec = 'double',
      compression = 6,
      chunksizes = c(dim(climatology$value)[1 : 2], 4)
    )

    log_debug('Calculating values')
    values <- lapply_variable(function(variable_index) {
      tensor::tensor(
        climatology$value[, , variable_index],
        X[, variable_index]
      )
    })
    residual <- input$value - Reduce(`+`, values)

    log_debug('Writing output to {output_filenames[file_index]}')
    output_fn <- ncdf4::nc_create(
      output_filenames[file_index],
      Filter(
        function(variable) {
          return (!has_fields || (variable$name %in% fields))
        },
        c(variables, list(residual_variable))
      ),
      force_v4 = TRUE
    )
    for (variable_index in seq_along(climatology$variable)) {
      variable_name <- climatology$variable[variable_index]
      if (has_fields && !(variable_name %in% fields)) {
        next
      }

      log_trace('Writing {climatology$variable[variable_index]}')
      ncdf4::ncvar_put(
        output_fn,
        variables[[variable_index]],
        values[[variable_index]]
      )
      ncdf4::nc_sync(output_fn)
    }
    log_trace('Writing residual')
    ncdf4::ncvar_put(
      output_fn,
      residual_variable,
      residual
    )
    log_trace('Closing')
    ncdf4::nc_close(output_fn)
  }
}

expand_climatology <- function(
  times,
  output_filename,
  climatology_filename,
  field_name,
  time_resolution = c('hourly', 'monthly'),
  time_alignment = c('start', 'middle', 'end')
) {
  time_resolution <- match.arg(time_resolution)
  time_alignment <- match.arg(time_alignment)

  log_debug('Loading climatology')
  climatology <- read_gridded_data(
    climatology_filename,
    'coefficient',
    include_time = FALSE,
    include_variable = TRUE
  )

  log_debug('Calculating design matrix over {length(times)} time periods')
  X <- climatology_design_matrix(
    time = times,
    origin = ymd_hms(attr(climatology, 'origin')),
    time_alignment = time_alignment,
    time_resolution = time_resolution,
    harmonics = attr(climatology, 'harmonics')
  )

  lapply_variable <- function(fn) {
    lapply(seq_along(climatology$variable), fn)
  }

  dimensions <- lapply(list(
    list(
      name = 'lon',
      units = 'degrees_east',
      longname = 'Longitude',
      vals = climatology$longitude
    ),
    list(
      name = 'lat',
      units = 'degrees_north',
      longname = 'Latitude',
      vals = climatology$latitude
    ),
    list(
      name = 'time',
      units = sprintf('hours since %s', attr(climatology, 'origin')),
      longname = 'Time',
      vals = as.integer(round(as.double(
        times - ymd_hms(attr(climatology, 'origin')),
        units = 'hours'
      )))
    )
  ), function(args) do.call(ncdf4::ncdim_def, args))

  log_trace('Defining variables')
  variables <- lapply_variable(function(variable_index) {
    variable_name <- climatology$variable[variable_index]
    ncdf4::ncvar_def(
      variable_name,
      units = climatology$value_units,
      longname = sprintf('%s %s', field_name, variable_name),
      dim = dimensions,
      prec = 'double',
      compression = 6,
      chunksizes = c(dim(climatology$value)[1 : 2], 4)
    )
  })

  log_debug('Calculating values')
  values <- lapply_variable(function(variable_index) {
    tensor::tensor(
      climatology$value[, , variable_index],
      X[, variable_index]
    )
  })

  log_debug('Writing output to {output_filename}')
  output_fn <- ncdf4::nc_create(output_filename, variables, force_v4 = TRUE)
  for (variable_index in seq_along(climatology$variable)) {
    log_trace('Writing {climatology$variable[variable_index]}')
    ncdf4::ncvar_put(
      output_fn,
      variables[[variable_index]],
      values[[variable_index]]
    )
    ncdf4::nc_sync(output_fn)
  }
  log_trace('Closing')
  ncdf4::nc_close(output_fn)
}


read_gridded_data <- function(
  filename,
  field_name,
  include_time = TRUE,
  include_variable = FALSE
) {
  with_nc_file(list(fn = filename), {
    output <- list(
      latitude = as.vector(ncdf4::ncvar_get(fn, 'lat')),
      longitude = as.vector(ncdf4::ncvar_get(fn, 'lon')),
      value = ncdf4::ncvar_get(fn, field_name),
      value_units = fn$var[[field_name]]$units,
      value_longname = fn$var[[field_name]]$longname
    )
    if (include_time) {
      output$time <- ncvar_get_time(fn, 'time')
    }
    if (include_variable) {
      output$variable <- as.vector(ncdf4::ncvar_get(fn, 'variable'))
    }
    attrs <- ncdf4::ncatt_get(fn, 0)
    for (name in names(attrs)) {
      attr(output, name) <- attrs[[name]]
    }
  })
  output
}

climatology_design_matrix <- function(
  time,
  origin,
  time_resolution = c('hourly', 'monthly'),
  time_alignment = c('start', 'middle', 'end'),
  harmonics
) {
  time_resolution <- match.arg(time_resolution)
  time_alignment <- match.arg(time_alignment)

  if (time_alignment != 'middle') {
    if (time_resolution == 'monthly') {
      time <- time + seconds(
        c('start' = 1, 'end' = -1)[time_alignment]
        * 24 * 60 * 60 * days_in_month(time) / 2
      )
    } else {
      time <- time + seconds(
        c('start' = 1, 'end' = -1)[time_alignment]
        * 60 * 30
      )
    }
  }

  trend <- as.numeric(
    time - origin,
    unit = 'hours'
  ) / (365.25 * 24)
  X <- matrix(NA, nrow = length(time), ncol = 4 * harmonics + 2)
  X[, 1] <- 1
  X[, 2] <- trend
  for (k in seq_len(harmonics)) {
    X[, 3 + 2 * (k - 1)] <- cospi(2 * trend / (1 / k))
    X[, 3 + 2 * (k - 1) + 1] <- sinpi(2 * trend / (1 / k))
    X[, 3 + 2 * harmonics + 2 * (k - 1)] <- trend * cospi(
      2 * trend / (1 / k)
    )
    X[, 3 + 2 * harmonics + 2 * (k - 1) + 1] <- trend * sinpi(
      2 * trend / (1 / k)
    )
  }
  colnames(X) <- c(
    'intercept',
    'trend',
    as.vector(sapply(seq_len(harmonics), function(k) {
      sprintf(c('cos12_%d_intercept', 'sin12_%d_intercept'), k)
    })),
    as.vector(sapply(seq_len(harmonics), function(k) {
      sprintf(c('cos12_%d_trend', 'sin12_%d_trend'), k)
    }))
  )
  X
}
