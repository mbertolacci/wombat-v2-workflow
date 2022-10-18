library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(ncdf4)

source(Sys.getenv('UTILS_PARTIAL'))

replace_na0 <- function(x) {
  x[is.na(x)] <- 0
  x
}

parser <- ArgumentParser()
parser$add_argument('--basis-vectors')
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations-augmented')
parser$add_argument('--samples')
parser$add_argument('--output')
args <- parser$parse_args()

basis_vectors <- fst::read_fst(args$basis_vectors)
control_emissions <- fst::read_fst(args$control_emissions) %>%
  filter(year(time) > 2014, year(time) < 2021)
perturbations <- fst::read_fst(args$perturbations_augmented) %>%
  filter(year(time) > 2014, year(time) < 2021) %>%
  mutate(
    minor_component = case_when(
      component == 'residual' ~ 'residual',
      component %in% c('intercept', 'trend') ~ 'linear',
      TRUE ~ 'periodic'
    )
  )
samples <- readRDS(args$samples)

alpha <- samples$alpha_df
n_samples <- 100L

times <- control_emissions %>%
  pull(time) %>%
  sort() %>%
  unique()
epoch <- min(times)
dimensions <- list(
  time = ncdim_def(
    'time',
    vals = as.numeric(times - epoch, units = 'days'),
    units = sprintf('days since %s', format(epoch, '%Y-%m-%d %H:%M:%S')),
    longname = 'time'
  ),
  latitude = ncdim_def(
    'latitude',
    vals = sort(unique(control_emissions$latitude)),
    units = 'degrees_north',
    longname = 'latitude'
  ),
  longitude = ncdim_def(
    'longitude',
    vals = sort(unique(control_emissions$longitude)),
    units = 'degrees_east',
    longname = 'longitude'
  ),
  sample = ncdim_def(
    'sample',
    vals = seq_len(n_samples),
    units = '',
    create_dimvar = FALSE
  )
)

output_grid <- control_emissions %>%
  distinct(longitude, latitude, time)

variable_longnames <- c(
  area = 'Area of grid cell',
  biofuel_total = 'Total flux for biofuels component',
  biomass_total = 'Total flux for biomass burning (i.e., wildfires) component',
  fossil_fuel_total = 'Total flux for fossil fuel component',
  respiration_linear = 'Linear component of respiration fluxes',
  respiration_periodic_bottom_up = 'Bottom-up estimate of time-varying periodic component of respiration fluxes',
  respiration_periodic_posterior = 'Posterior samples for time-varying periodic component of respiration fluxes',
  respiration_residual_bottom_up = 'Bottom-up estimate of time-varying residual component of respiration fluxes',
  respiration_residual_posterior = 'Posterior samples for time-varying residual component of respiration fluxes',
  gpp_linear_bottom_up = 'Bottom-up estimate of time-varying linear component of GPP fluxes',
  gpp_linear_posterior = 'Posterior samples for time-varying linear component of GPP fluxes',
  gpp_periodic_bottom_up = 'Bottom-up estimate of time-varying periodic component of GPP fluxes',
  gpp_periodic_posterior = 'Posterior samples for time-varying periodic component of GPP fluxes',
  gpp_residual_bottom_up = 'Bottom-up estimate of time-varying residual component of GPP fluxes',
  gpp_residual_posterior = 'Posterior samples for time-varying residual component of GPP fluxes',
  ocean_linear = 'Linear component of ocean fluxes',
  ocean_periodic = 'Time-varying periodic component of ocean fluxes',
  ocean_residual_bottom_up = 'Bottom-up estimate of time-varying residual component of ocean fluxes',
  ocean_residual_posterior = 'Posterior samples for time-varying residual component of ocean fluxes'
)

variables <- list()
variables[['area']] <- ncvar_def(
  'area',
  units = 'm2',
  dim = dimensions[c('longitude', 'latitude')]
)
# Known values
for (name in c(
  'respiration_linear',
  'ocean_linear',
  'ocean_periodic',
  'biofuel_total',
  'biomass_total',
  'fossil_fuel_total',
  'respiration_periodic_bottom_up',
  'respiration_residual_bottom_up',
  'gpp_linear_bottom_up',
  'gpp_periodic_bottom_up',
  'gpp_residual_bottom_up',
  'ocean_residual_bottom_up'
)) {
  variables[[name]] <- ncvar_def(
    name,
    units = 'kgCO2/m2/s',
    dim = dimensions[c('longitude', 'latitude', 'time')],
    compression = 9
  )
}
# Unknown values
for (name in c(
  'respiration_periodic_posterior',
  'respiration_residual_posterior',
  'gpp_linear_posterior',
  'gpp_periodic_posterior',
  'gpp_residual_posterior',
  'ocean_residual_posterior'
)) {
  variables[[name]] <- ncvar_def(
    name,
    units = 'kgCO2/m2/s',
    dim = dimensions[c('sample', 'longitude', 'latitude', 'time')],
    compression = 9
  )
}

for (name_i in names(variables)) {
  variables[[name_i]]$longname <- variable_longnames[name_i]
}

output_name_to_input <- c(
  'respiration' = 'bio_resp_tot',
  'gpp' = 'bio_assim',
  'ocean' = 'ocean'
)

output_fn <- nc_create(args$output, variables, force_v4 = TRUE)
ncvar_put(
  output_fn,
  'area',
  control_emissions %>%
    distinct(longitude, latitude, area) %>%
    pull(area) %>%
    array(dim = variables[['area']]$varsize)
)
# Known values in the control emissions
for (inventory_i in c('biofuel', 'biomass', 'fossil_fuel')) {
  name_i <- sprintf('%s_total', inventory_i)
  log_trace('Outputting {name_i}')
  ncvar_put(
    output_fn,
    name_i,
    control_emissions %>%
      filter(inventory == inventory_i) %>%
      pull(value) %>%
      array(dim = variables[[name_i]]$varsize)
  )
}
# Known values in the perturbations
for (name_i in c(
  'respiration_linear',
  'ocean_linear',
  'ocean_periodic',
  'respiration_periodic_bottom_up',
  'respiration_residual_bottom_up',
  'gpp_linear_bottom_up',
  'gpp_periodic_bottom_up',
  'gpp_residual_bottom_up',
  'ocean_residual_bottom_up'
)) {
  parts_i <- strsplit(name_i, '_')[[1]]
  inventory_i <- output_name_to_input[parts_i[1]]
  minor_component_i <- parts_i[2]

  log_trace('Outputting {name_i} {inventory_i} {minor_component_i}')

  output_values_i <- output_grid %>%
    left_join(
      perturbations %>%
        filter(inventory == inventory_i, minor_component == minor_component_i) %>%
        group_by(longitude, latitude, time) %>%
        summarise(value = sum(value)),
      by = c('longitude', 'latitude', 'time')
    ) %>%
    mutate(value = replace_na0(value))

  ncvar_put(
    output_fn,
    name_i,
    array(output_values_i$value, dim = variables[[name_i]]$varsize)
  )
}

# Unknown values
for (name_i in c(
  'respiration_periodic_posterior',
  'respiration_residual_posterior',
  'gpp_linear_posterior',
  'gpp_periodic_posterior',
  'gpp_residual_posterior',
  'ocean_residual_posterior'
)) {
  parts_i <- strsplit(name_i, '_')[[1]]
  inventory_i <- output_name_to_input[parts_i[1]]
  minor_component_i <- parts_i[2]

  log_trace('Outputting {name_i} {inventory_i} {minor_component_i}')

  samples_i <- perturbations %>%
    filter(inventory == inventory_i, minor_component == minor_component_i) %>%
    left_join(
      alpha %>%
        select(basis_vector, alpha = value, alpha_samples = value_samples) %>%
        mutate(
          alpha_samples = alpha_samples[, floor(seq(1L, 500L, length.out = n_samples))]
        ),
      by = 'basis_vector'
    ) %>%
    rename(perturbation = value) %>%
    mutate(
      value_samples = perturbation * (1 + replace_na0(alpha_samples))
    ) %>%
    select(-perturbation) %>%
    group_by(longitude, latitude, time) %>%
    summarise(
      value_samples = t(colSums(value_samples))
    )

  output_values_i <- output_grid %>%
    left_join(
      samples_i,
      by = c('longitude', 'latitude', 'time')
    ) %>%
    mutate(value_samples = replace_na0(value_samples))

  ncvar_put(
    output_fn,
    name_i,
    array(as.vector(t(output_values_i$value_samples)), dim = variables[[name_i]]$varsize)
  )
}
ncatt_put(output_fn, 0, 'contact', 'Michael Bertolacci <m.bertolacci@gmail.com>')
ncatt_put(output_fn, 0, 'creation_date', format(today()))
nc_close(output_fn)
