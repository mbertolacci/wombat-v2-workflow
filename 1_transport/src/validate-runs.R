library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(readr)
library(lubridate, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--runs')
args <- parser$parse_args()

list_run_outputs <- function(paths, pattern) {
  sort(do.call(c, lapply(paths, function(path) {
    list.files(
      file.path(path, 'output'),
      pattern = pattern,
      full.names = TRUE
    )
  })))
}

list_species_conc_files <- function(paths, resolution) {
  list_run_outputs(paths, sprintf(
    'GEOSChem\\.SpeciesConc%s\\..*',
    paste0(toupper(substring(resolution, 1, 1)), substring(resolution, 2))
  ))
}

run_check <- function(run, lhs, rhs) {
  lhs_q <- deparse(substitute(lhs))
  rhs_q <- deparse(substitute(rhs))
  if (lhs == rhs) return(FALSE)
  log_warn('[{run}] {lhs_q} == {rhs_q} not satisfied; {lhs} != {rhs}')
  return(TRUE)
}

mapping <- read_csv(file.path(args$runs, 'mapping.csv'), show_col_types = FALSE)
all_runs <- list.files(args$runs)

log_info('Validating runs')
any_issue <- FALSE
for (run in c('base', sort(unique(mapping$run)))) {
  run_parts <- sort(file.path(
    args$runs,
    all_runs[startsWith(all_runs, run)]
  ))
  run_start_time <- NULL
  run_end_time <- NULL
  for (run_part in run_parts) {
    config_lines <- readLines(file.path(run_part, 'input.geos'), n = 10)
    run_start_time_current <- ymd_hms(config_lines[startsWith(config_lines, 'Start')])
    run_end_time_current <- ymd_hms(config_lines[startsWith(config_lines, 'End')])
    if (is.null(run_start_time) || run_start_time > run_start_time_current) {
      run_start_time <- run_start_time_current
    }
    if (is.null(run_end_time) || run_end_time < run_end_time_current) {
      run_end_time <- run_end_time_current
    }
  }

  n_days <- interval(run_start_time, run_end_time) %/% days(1)
  n_initial_days <- interval(
    run_start_time,
    min(run_end_time, run_start_time + months(4))
  ) %/% days(1)
  n_months <- interval(run_start_time, run_end_time) %/% months(1)

  hourly_species_conc_files <- list_species_conc_files(run_parts, 'hourly')
  daily_species_conc_files <- list_species_conc_files(run_parts, 'daily')
  hemco_diagnostics <- list_run_outputs(run_parts, 'HEMCO_diagnostics.*')

  any_issue <- run_check(run, length(daily_species_conc_files), n_months) || any_issue
  any_issue <- run_check(run, length(hemco_diagnostics), n_months) || any_issue
  if (startsWith(run, 'residual')) {
    any_issue <- run_check(run, length(hourly_species_conc_files), n_initial_days) || any_issue
  } else {
    any_issue <- run_check(run, length(hourly_species_conc_files), n_days) || any_issue
  }
}

if (any_issue) {
  stop('Run validation failed')
}

log_info('Done')
