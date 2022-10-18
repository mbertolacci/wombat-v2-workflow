library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(fst)
library(parallel)

source(Sys.getenv('UTILS_PARTIAL'))

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

mclapply_chunked <- function(X, FUN, chunk_size, ...){
  chunk_size <- if (is.numeric(chunk_size) && length(chunk_size) == 1) {
    as.integer(chunk_size)
  } else {
    1L
  }
  if (chunk_size == 1L) return(lapply(X, FUN, ...))
  N <- length(X)
  index_list <- splitIndices(N, N / chunk_size)
  result_list <- list()
  for (i in seq_along(index_list)){
    log_trace('[mclapply_chunked] Processing chunk {i} / {length(index_list)}')
    index_vec <- index_list[[i]]
    gc()

    parts <- mclapply(X[index_vec], FUN, ...)
    for (part in parts) {
      if (is(part, 'try-error')) {
        print(part)
        stop('Had a failure')
      }
    }
    result_list[index_vec] <- parts
    rm(parts)
    gc()
  }
  result_list
}

parser <- ArgumentParser()
parser$add_argument('--input', nargs = '+')
parser$add_argument('--resolution', nargs = '+')
parser$add_argument('--runs')
parser$add_argument('--matched-runs')
parser$add_argument('--output-base', nargs = '+')
args <- parser$parse_args()

mapping <- readr::read_csv(
  file.path(args$runs, 'mapping.csv'),
  show_col_types = FALSE
)
mapping_background <- mapping %>%
  filter(startsWith(basis_function, 'background')) %>%
  mutate(
    run_base = sapply(strsplit(run, '_part'), getElement, 1)
  )
mapping_basis_function <- mapping %>%
  filter(!startsWith(basis_function, 'background')) %>%
  mutate(
    species = factor(species),
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

inputs <- do.call(
  vctrs::vec_recycle_common,
  args[c('input', 'resolution', 'output_base')]
)
runs <- unique(mapping$run)

for (index in seq_along(inputs$input)) {
  log_info('Loading parts for {inputs$input[index]} and computing sensitivities')

  log_debug('Loading background values')
  background_values <- mclapply(mapping_background$run, function(run_i) {
    log_trace('Loading background from {run_i}')
    read_fst(
      file.path(args$matched_runs, run_i, inputs$input[index]),
      columns = c('species', 'value')
    ) %>%
      filter(species %in% mapping_background$species) %>%
      pull(value)
  }, mc.cores = get_cores())
  names(background_values) <- mapping_background$run_base

  log_debug('Counting values per run')
  n_per_run <- simplify2array(mclapply(runs, function(run_i) {
    read_fst(
      file.path(args$matched_runs, run_i, inputs$input[index]),
      columns = 'species'
    ) %>%
      filter(!(species %in% mapping_background$species)) %>%
      nrow()
  }, mc.cores = get_cores()))

  n_max <- 2^31 - 1
  part_indices <- NULL
  n_current <- 0
  indices_current <- NULL
  for (i in seq_along(n_per_run)) {
    if (n_current + n_per_run[i] > n_max) {
      part_indices <- c(part_indices, list(indices_current))
      n_current <- 0
      indices_current <- NULL
    }
    n_current <- n_current + n_per_run[i]
    indices_current <- c(indices_current, i)
  }
  part_indices <- c(part_indices, list(indices_current))

  for (part_index in seq_along(part_indices)) {
    log_debug('Loading part {part_index} / {length(part_indices)}')
    output_parts <- mclapply_chunked(runs[part_indices[[part_index]]], function(run_i) {
      log_debug('Processing {run_i}')

      filename <- file.path(args$matched_runs, run_i, inputs$input[index])
      log_trace('Reading {filename}')
      tracers <- read_fst(
        filename,
        columns = c('observation_id', 'species', 'value')
      )
      mapping_basis_function_i <- mapping_basis_function %>%
        filter(run == run_i)
      log_trace('Computing sensitivities')
      tracers %>%
        filter(
          species %in% mapping_basis_function_i$species
        ) %>%
        group_by(species) %>%
        mutate(
          value = value - background_values[[
            strsplit(run_i, '_part')[[1]][1]
          ]]
        ) %>%
        ungroup() %>%
        filter(value != 0) %>%
        left_join(
          mapping_basis_function_i %>%
            select(-run),
          by = 'species'
        ) %>%
        mutate(
          resolution = factor(
            inputs$resolution[index],
            levels = c('daily', 'hourly')
          )
        ) %>%
        select(observation_id, resolution, region, month, inventory, component, value)
    }, mc.cores = get_cores(), chunk_size = 2 * get_cores())

    log_debug('Combining')
    output <- bind_rows(output_parts)

    output_filename <- sprintf('%s-part-%d.fst', inputs$output_base[index], part_index)
    log_info('Writing to {output_filename}')
    fst::write_fst(output, output_filename)

    log_info('Done')
    rm(output)
    rm(output_parts)
    print(gc())
  }
}
