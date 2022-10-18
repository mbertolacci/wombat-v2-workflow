library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(Matrix)
library(ncdf4)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--control-emissions')
parser$add_argument('--region-grid')
parser$add_argument('--samples')
parser$add_argument('--sib4-climatology-assim')
parser$add_argument('--sib4-climatology-resp-tot')
parser$add_argument('--landschutzer-climatology')
parser$add_argument('--output')
args <- parser$parse_args()

log_debug('Loading inputs')
control_emissions <- fst::read_fst(args$control_emissions)
samples <- readRDS(args$samples)
region_grid <- readRDS(args$region_grid)

fit_region_grid_df <- NULL
for (latitude_index in seq_along(attr(region_grid$Region01, 'latitude'))) {
  for (longitude_index in seq_along(attr(region_grid$Region01, 'longitude'))) {
    found_any <- FALSE
    for (region_name in names(region_grid)) {
      if (region_grid[[region_name]][longitude_index, latitude_index] > 0) {
        fit_region_grid_df <- rbind(fit_region_grid_df, data.frame(
          longitude = attr(region_grid$Region01, 'longitude')[longitude_index],
          latitude = attr(region_grid$Region01, 'latitude')[latitude_index],
          region = region_name
        ))
        found_any <- TRUE
      }
    }
    if (!found_any) {
      fit_region_grid_df <- rbind(fit_region_grid_df, data.frame(
        longitude = attr(region_grid$Region01, 'longitude')[longitude_index],
        latitude = attr(region_grid$Region01, 'latitude')[latitude_index],
        region = 'Region00'
      ))
    }
  }
}

X_fit_region_to_2x25 <- with(
  fit_region_grid_df %>%
    mutate(region = factor(region)),
  sparseMatrix(
    i = seq_along(longitude),
    j = as.integer(region),
    x = rep(1, length(longitude)),
    dims = c(length(longitude), nlevels(region))
  )
)
# NOTE(mgnb): first region is Region00, the nothing region
X_fit_region_to_2x25 <- X_fit_region_to_2x25[, 2 : ncol(X_fit_region_to_2x25)]

climatologies <- lapply(
  args[c('sib4_climatology_assim', 'sib4_climatology_resp_tot', 'landschutzer_climatology')],
  read_gridded_data,
  'coefficient',
  include_time = FALSE,
  include_variable = TRUE
)
origins <- sapply(climatologies, attr, 'origin')
harmonics <- sapply(climatologies, attr, 'harmonics')
stopifnot(all(origins == origins[1]))
origin <- origins[1]

climatology_prior <- bind_rows(
  read_climatology(args$sib4_climatology_assim) %>%
    mutate(
      value = -value,
      inventory = 'bio_assim'
    ),
  read_climatology(args$sib4_climatology_resp_tot) %>%
    mutate(inventory = 'bio_resp_tot'),
  read_climatology(args$landschutzer_climatology) %>%
    mutate(inventory = 'ocean')
)

log_debug('Computing posterior functionals')
n_samples <- 100
alpha_modified <- samples$alpha_df %>%
  select(
    inventory,
    region,
    variable = component,
    alpha_mean = value,
    alpha_samples = value_samples
  ) %>%
  mutate(
    variable = factor(as.character(variable), levels(climatology_prior$variable)),
    alpha_samples = alpha_samples[, floor(seq(1, ncol(alpha_samples), length.out = n_samples))]
  ) %>%
  filter(variable != 'residual')

climatology_samples <- climatology_prior %>%
  group_by(variable, inventory) %>%
  group_modify(~ {
    alpha_modified_i <- alpha_modified %>%
      filter(variable == .y$variable, inventory == .y$inventory)

    if (nrow(alpha_modified_i) == 0) {
      .x %>%
        mutate(
          value_prior = value,
          value = 0,
          value_samples = matrix(0, nrow = n(), ncol = n_samples)
        )
    } else {
      X_i <- Diagonal(x = .x$value) %*% X_fit_region_to_2x25[
        ,
        as.integer(alpha_modified_i$region)
      ]

      .x %>%
        mutate(
          value_prior = value,
          value_samples = value_prior + as.matrix(X_i %*% alpha_modified_i$alpha_samples),
          value = rowMeans(value_samples)
        )
    }
  }) %>%
  mutate(
    inventory = factor(
      c(
        'bio_assim' = 'GPP',
        'bio_resp_tot' = 'Respiration',
        'ocean' = 'Ocean'
      )[inventory],
      levels = c('GPP', 'Respiration', 'Ocean')
    )
  ) %>%
  arrange(latitude, longitude, inventory, variable)

dimensions <- list(
  latitude = ncdim_def(
    'latitude',
    vals = sort(unique(control_emissions$latitude)),
    units = 'degrees_north',
    longname = 'Latitude'
  ),
  longitude = ncdim_def(
    'longitude',
    vals = sort(unique(control_emissions$longitude)),
    units = 'degrees_east',
    longname = 'Longitude'
  ),
  inventory = ncdim_def(
    'inventory',
    vals = seq_along(levels(climatology_samples$inventory)),
    longname = 'Inventory',
    units = '',
    create_dimvar = FALSE
  ),
  sample = ncdim_def(
    'sample',
    vals = seq_len(n_samples),
    longname = 'Sample',
    units = '',
    create_dimvar = FALSE
  )
)

variable_longnames <- c(
  area = 'Area of grid cell',
  inventory = 'Name of inventory',
  intercept = 'flux intercept',
  trend = 'flux trend',
  cos12_1_intercept = 'intercept of 12-month cosine harmonic of flux',
  sin12_1_intercept = 'intercept of 12-month sine harmonic of flux',
  cos12_2_intercept = 'intercept of 6-month (12 / 2) cosine harmonic of flux',
  sin12_2_intercept = 'intercept of 6-month (12 / 2) sine harmonic of flux',
  cos12_3_intercept = 'intercept of 4-month (12 / 3) cosine harmonic of flux',
  sin12_3_intercept = 'intercept of 4-month (12 / 3) sine harmonic of flux',
  cos12_1_trend = 'trend of 12-month cosine harmonic of flux',
  sin12_1_trend = 'trend of 12-month sine harmonic of flux',
  cos12_2_trend = 'trend of 6-month (12 / 2) cosine harmonic of flux',
  sin12_2_trend = 'trend of 6-month (12 / 2) sine harmonic of flux',
  cos12_3_trend = 'trend of 4-month (12 / 3) cosine harmonic of flux',
  sin12_3_trend = 'trend of 4-month (12 / 3) sine harmonic of flux'
)

variables <- list()
variables[['area']] <- ncvar_def(
  'area',
  units = 'm2',
  longname = variable_longnames[['area']],
  dim = dimensions[c('longitude', 'latitude')]
)
variables[['inventory']] <- ncdf4::ncvar_def(
  'inventory',
  units = '',
  longname = variable_longnames[['inventory']],
  dim = list(
    ncdf4::ncdim_def(
      'nchar_inventory',
      '',
      seq_len(max(nchar(levels(climatology_samples$inventory)))),
      create_dimvar = FALSE
    ),
    dimensions[['inventory']]
  ),
  prec = 'char',
  compression = 9
)
for (name in c(
  'intercept',
  'trend',
  'cos12_1_intercept',
  'sin12_1_intercept',
  'cos12_2_intercept',
  'sin12_2_intercept',
  'cos12_3_intercept',
  'sin12_3_intercept',
  'cos12_1_trend',
  'sin12_1_trend',
  'cos12_2_trend',
  'sin12_2_trend',
  'cos12_3_trend',
  'sin12_3_trend'
)) {
  unit <- if (endsWith(name, 'trend')) 'kgCO2/m2/s/year' else 'kgCO2/m2/s'
  variables[[sprintf('%s_bottom_up', name)]] <- ncvar_def(
    sprintf('%s_bottom_up', name),
    units = unit,
    longname = sprintf('Bottom-up estimate of %s', variable_longnames[[name]]),
    dim = dimensions[c('longitude', 'latitude', 'inventory')],
    compression = 9
  )
  variables[[sprintf('%s_posterior_samples', name)]] <- ncvar_def(
    sprintf('%s_posterior_samples', name),
    units = unit,
    longname = sprintf('Posterior samples of %s', variable_longnames[[name]]),
    dim = dimensions[c('longitude', 'latitude', 'inventory', 'sample')],
    compression = 9
  )
}

output_fn <- nc_create(
  args$output,
  variables,
  force_v4 = TRUE
)
ncvar_put(
  output_fn,
  'area',
  control_emissions %>%
    distinct(longitude, latitude, area) %>%
    pull(area) %>%
    array(dim = variables[['area']]$varsize)
)
value_dim_base <- c(
  length(unique(control_emissions$longitude)),
  length(unique(control_emissions$latitude)),
  nlevels(climatology_samples$inventory)
)
for (variable_i in levels(climatology_samples$variable)) {
  log_trace('Writing {variable_i}')
  value_bottom_up_i <- array(0, dim = value_dim_base)
  value_posterior_samples_i <- array(0, dim = c(value_dim_base, n_samples))

  for (j in seq_len(nlevels(climatology_samples$inventory))) {
    climatology_samples_i_j <- climatology_samples %>%
      filter(
        variable == variable_i,
        inventory == levels(inventory)[j]
      )
    if (nrow(climatology_samples_i_j) == 0) next

    value_bottom_up_i[, , j] <- climatology_samples_i_j$value_prior
    value_posterior_samples_i[, , j, ] <- climatology_samples_i_j$value_samples
  }
  ncvar_put(
    output_fn,
    sprintf('%s_bottom_up', variable_i),
    value_bottom_up_i
  )
  ncvar_put(
    output_fn,
    sprintf('%s_posterior_samples', variable_i),
    value_posterior_samples_i
  )
}
ncvar_put(output_fn, 'inventory', levels(climatology_samples$inventory))
ncatt_put(output_fn, 0, 'origin', origin)
ncatt_put(output_fn, 0, 'contact', 'Michael Bertolacci <m.bertolacci@gmail.com>')
ncatt_put(output_fn, 0, 'creation_date', format(today()))
nc_close(output_fn)
