library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(ncdf4)
library(parallel)
library(lubridate, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))
Rcpp::sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))

strip_to_starts <- function(x, starts) {
  output <- character(length(x))
  for (start in starts) {
    output[startsWith(x, start)] <- start
  }
  output
}

strip_off_starts <- function(x, starts) {
  output <- character(length(x))
  for (start in starts) {
    which_start <- startsWith(x, start)
    output[which_start] <- substring(
      x[which_start],
      nchar(start) + 1
    )
  }
  output
}

parser <- ArgumentParser()
parser$add_argument('--runs', nargs = '+')
parser$add_argument('--matched-runs', nargs = '+')
parser$add_argument('--output')
args <- parser$parse_args()

stopifnot(length(args$runs) == length(args$matched_runs))

outputs <- lapply(seq_along(args$runs), function(index) {
  log_debug('Reading mapping {index}')
  mapping_basis_function <- readr::read_csv(
    file.path(args$runs[index], 'mapping.csv'),
    show_col_types = FALSE
  ) %>%
    filter(!startsWith(basis_function, 'background')) %>%
    mutate(
      parts = strsplit(basis_function, '_region'),
      part1 = sapply(parts, getElement, 1),
      parts2 = strsplit(sapply(parts, getElement, 2), '_month'),
      region = factor(sapply(parts2, getElement, 1)),
      month = factor(sapply(parts2, function(x) {
        if (length(x) == 2) x[[2]] else NA
      })),
      inventory = factor(strip_to_starts(
        part1,
        c('ocean', 'bio_assim', 'bio_resp_tot')
      )),
      component = factor(strip_off_starts(
        part1,
        c('ocean_', 'bio_assim_', 'bio_resp_tot_')
      ))
    ) %>%
    select(-basis_function, -parts, -part1, -parts2)

  runs <- unique(mapping_basis_function$run)
  mclapply(runs, function(run_i) {
    log_debug('Processing {run_i}')

    filename <- file.path(args$matched_runs[index], run_i, 'monthly-fluxes.nc')
    fn <- nc_open(filename)
    on.exit(nc_close(fn))
    v <- function(...) ncvar_get(fn, ...)

    locations <- expand.grid(
      longitude = as.vector(v('lon')),
      latitude = as.vector(v('lat')),
      time = as.Date(ncvar_get_time(fn, 'time'))
    )
    mapping_i <- mapping_basis_function %>% filter(run == run_i)
    bind_rows(lapply(mapping_i$species, function(species_i) {
      as_tibble(cbind(locations, data.frame(
        species = species_i,
        value = as.vector(v(sprintf('Emis_%s', species_i)))
      )))
    })) %>%
      left_join(
        mapping_i %>%
          select(species, region, month, inventory, component),
        by = 'species'
      ) %>%
      select(-species)
  }, mc.cores = get_cores()) %>%
    bind_rows()
})

log_debug('Merging parts')
perturbations <- outputs[[1]]
if (length(args$runs) > 1) {
  for (index in 2 : length(args$runs)) {
    perturbations <- perturbations %>%
      left_join(
        outputs[[index]] %>%
          mutate(
            region = update_levels(region, levels(perturbations$region)),
            month = update_levels(month, levels(perturbations$month)),
            inventory = update_levels(inventory, levels(perturbations$inventory)),
            component = update_levels(component, levels(perturbations$component))
          ) %>%
          rename(value2 = value),
        by = c(
          'region',
          'month',
          'inventory',
          'component',
          'latitude',
          'longitude',
          'time'
        )
      ) %>%
      mutate(
        value = ifelse(is.na(value2), value, value2)
      ) %>%
      select(-value2)
  }
}

log_debug('Removing zero values')
perturbations <- perturbations %>%
  filter(value != 0)

log_debug('Writing to {args$output}')
fst::write_fst(perturbations, args$output)

log_debug('Done')
