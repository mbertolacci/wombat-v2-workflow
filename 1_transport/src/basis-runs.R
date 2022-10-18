library(argparse)
library(lubridate, warn.conflicts = FALSE)
library(wombatbasis)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--template-directory')
parser$add_argument('--region-mask')
parser$add_argument('--region', nargs = '+')
parser$add_argument('--fossil-inventory-path')
parser$add_argument('--fires-inventory-path')
parser$add_argument('--ocean-climatology-inventory-path')
parser$add_argument('--ocean-residual-path')
parser$add_argument('--sib4-assim-climatology-inventory-path')
parser$add_argument('--sib4-assim-residual-path')
parser$add_argument('--sib4-resp-tot-climatology-inventory-path')
parser$add_argument('--sib4-resp-tot-residual-path')
parser$add_argument('--geos-chem-restart')
parser$add_argument('--geos-chem-extdata')
parser$add_argument('--geos-chem-source')
parser$add_argument('--output-directory')
args <- parser$parse_args()

source(file.path(args$template_directory, 'base-run.R'))

offset_to_root <- function(x) file.path(root_offset, x)

root_offset <- '../../../..'
base_run <- build_base_run(
  offset_to_root(args$geos_chem_extdata),
  offset_to_root(args$fossil_inventory_path),
  offset_to_root(args$fires_inventory_path),
  offset_to_root(args$ocean_climatology_inventory_path),
  offset_to_root(args$ocean_residual_path),
  offset_to_root(args$sib4_assim_climatology_inventory_path),
  offset_to_root(args$sib4_assim_residual_path),
  offset_to_root(args$sib4_resp_tot_climatology_inventory_path),
  offset_to_root(args$sib4_resp_tot_residual_path)
)
base_run$symlinks[['CodeDir']] <- offset_to_root(args$geos_chem_source)

log_info('Reading regions')
region_grid <- readRDS(args$region_mask)
if (length(args[['region']]) > 0) {
  region_grid <- region_grid[args[['region']]]
}

months <- seq(as.Date('2014-09-01'), as.Date('2021-03-01'), by = 'month')
names(months) <- format(months, '%Y-%m')
month_factor_grid <- times_to_grids(months)

log_info('Constructing basis')
basis <- flux_basis(
  structure = ~ (
    (OCEAN_LSCHULZ_INTERCEPT | region | ocean_intercept)
    + (OCEAN_LSCHULZ_TREND | region | ocean_trend)
    + (OCEAN_LSCHULZ_COS12_1_INTERCEPT | region | ocean_cos12_1_intercept)
    + (OCEAN_LSCHULZ_SIN12_1_INTERCEPT | region | ocean_sin12_1_intercept)
    + (OCEAN_LSCHULZ_COS12_2_INTERCEPT | region | ocean_cos12_2_intercept)
    + (OCEAN_LSCHULZ_SIN12_2_INTERCEPT | region | ocean_sin12_2_intercept)
    + (OCEAN_LSCHULZ_COS12_1_TREND | region | ocean_cos12_1_trend)
    + (OCEAN_LSCHULZ_SIN12_1_TREND | region | ocean_sin12_1_trend)
    + (OCEAN_LSCHULZ_COS12_2_TREND | region | ocean_cos12_2_trend)
    + (OCEAN_LSCHULZ_SIN12_2_TREND | region | ocean_sin12_2_trend)
    + (SIB4_ASSIM_INTERCEPT | region | bio_assim_intercept)
    + (SIB4_ASSIM_TREND | region | bio_assim_trend)
    + (SIB4_ASSIM_COS12_1_INTERCEPT | region | bio_assim_cos12_1_intercept)
    + (SIB4_ASSIM_SIN12_1_INTERCEPT | region | bio_assim_sin12_1_intercept)
    + (SIB4_ASSIM_COS12_2_INTERCEPT | region | bio_assim_cos12_2_intercept)
    + (SIB4_ASSIM_SIN12_2_INTERCEPT | region | bio_assim_sin12_2_intercept)
    + (SIB4_ASSIM_COS12_3_INTERCEPT | region | bio_assim_cos12_3_intercept)
    + (SIB4_ASSIM_SIN12_3_INTERCEPT | region | bio_assim_sin12_3_intercept)
    + (SIB4_ASSIM_COS12_1_TREND | region | bio_assim_cos12_1_trend)
    + (SIB4_ASSIM_SIN12_1_TREND | region | bio_assim_sin12_1_trend)
    + (SIB4_ASSIM_COS12_2_TREND | region | bio_assim_cos12_2_trend)
    + (SIB4_ASSIM_SIN12_2_TREND | region | bio_assim_sin12_2_trend)
    + (SIB4_ASSIM_COS12_3_TREND | region | bio_assim_cos12_3_trend)
    + (SIB4_ASSIM_SIN12_3_TREND | region | bio_assim_sin12_3_trend)
    + (SIB4_RESP_TOT_INTERCEPT | region | bio_resp_tot_intercept)
    + (SIB4_RESP_TOT_TREND | region | bio_resp_tot_trend)
    + (SIB4_RESP_TOT_COS12_1_INTERCEPT | region | bio_resp_tot_cos12_1_intercept)
    + (SIB4_RESP_TOT_SIN12_1_INTERCEPT | region | bio_resp_tot_sin12_1_intercept)
    + (SIB4_RESP_TOT_COS12_2_INTERCEPT | region | bio_resp_tot_cos12_2_intercept)
    + (SIB4_RESP_TOT_SIN12_2_INTERCEPT | region | bio_resp_tot_sin12_2_intercept)
    + (SIB4_RESP_TOT_COS12_3_INTERCEPT | region | bio_resp_tot_cos12_3_intercept)
    + (SIB4_RESP_TOT_SIN12_3_INTERCEPT | region | bio_resp_tot_sin12_3_intercept)
    + (SIB4_RESP_TOT_COS12_1_TREND | region | bio_resp_tot_cos12_1_trend)
    + (SIB4_RESP_TOT_SIN12_1_TREND | region | bio_resp_tot_sin12_1_trend)
    + (SIB4_RESP_TOT_COS12_2_TREND | region | bio_resp_tot_cos12_2_trend)
    + (SIB4_RESP_TOT_SIN12_2_TREND | region | bio_resp_tot_sin12_2_trend)
    + (SIB4_RESP_TOT_COS12_3_TREND | region | bio_resp_tot_cos12_3_trend)
    + (SIB4_RESP_TOT_SIN12_3_TREND | region | bio_resp_tot_sin12_3_trend)
    + (OCEAN_LSCHULZ_RESIDUAL | region * month | ocean_residual)
    + (SIB4_ASSIM_RESIDUAL | region * month | bio_assim_residual)
    + (SIB4_RESP_TOT_RESIDUAL | region * month | bio_resp_tot_residual)
  ),
  base_run = base_run,
  factor_grids = list(
    region = region_grid,
    month = month_factor_grid
  )
)

read_template_file <- function(filename) {
  paste0(
    c(readLines(file.path(args$template_directory, filename)), ''),
    collapse = '\n'
  )
}

postprocess_script <- read_template_file('postprocess-run.sh')

log_info('Converting basis to runs')
runs <- basis_to_runs(
  basis,
  groups = sapply(basis$basis_functions, function(basis_function) {
    if (endsWith(basis_function$part, 'residual')) {
      'residual'
    } else {
      'climatology'
    }
  }),
  strategy = 'passive_species',
  max_species_per_run = 24,
  max_run_length = years(1),
  postprocess_run = function(run) {
    if (run$name == 'base') {
      run$configuration$symlinks[[
        basename(args$geos_chem_restart)
      ]] <- offset_to_root(args$geos_chem_restart)
    }

    start_date <- run$configuration$main$simulation_menu$start
    end_date <- run$configuration$main$simulation_menu$end

    hourly_start_to_remove <- start_date + months(4)
    if (
      run$name == 'base'
      || startsWith(run$name, 'climatology')
      || hourly_start_to_remove > end_date
    ) {
      hourly_months_to_remove <- ''
    } else {
      hourly_months_to_remove <- paste0(format(seq.Date(
        hourly_start_to_remove,
        end_date,
        by = 'month'
      ), '%Y%m'), collapse = ' ')
    }

    hemco_start_to_remove <- start_date + months(1)
    if (
      run$name == 'base'
      || startsWith(run$name, 'climatology')
      || hemco_start_to_remove > end_date
    ) {
      hemco_months_to_remove <- ''
    } else {
      hemco_months_to_remove <- paste0(format(seq.Date(
        hemco_start_to_remove,
        end_date,
        by = 'month'
      ), '%Y%m'), collapse = ' ')
    }

    if (run$name != 'base') {
      run$configuration$history$collections <- Filter(function(collection) {
        return(!(
          startsWith(collection$name, 'LevelEdge')
          || startsWith(collection$name, 'StateMet')
        ))
      }, run$configuration$history$collections)
    }

    run$configuration$files[['postprocess-run.sh']] <- sprintf(
      postprocess_script,
      hemco_months_to_remove,
      if (run$name == 'base') {
        'LevelEdge StateMet SpeciesConc'
      } else {
        'SpeciesConc'
      },
      hourly_months_to_remove
    )
    for (filename in c('run-niasra-hpc-sbatch.sh', 'run-gadi.sh')) {
      run$configuration$files[[filename]] <- read_template_file(filename)
    }
    run
  }
)

log_info('Writing runs')
write_basis_runs(runs, args$output_directory)

log_info('Done')
