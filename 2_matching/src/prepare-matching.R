library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(readr)
library(glue, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--oco2-soundings')
parser$add_argument('--tccon-sounding-directory')
parser$add_argument('--obspack-directory')
parser$add_argument('--lauder')
parser$add_argument('--input')
parser$add_argument('--template-directory')
parser$add_argument('--output')
args <- parser$parse_args()

mapping <- read_csv(
  file.path(args$input, 'mapping.csv'),
  show_col_types = FALSE
)

read_template_file <- function(filename, data) {
  glue_data(
    data,
    paste0(
      c(readLines(file.path(args$template_directory, filename)), ''),
      collapse = '\n'
    ),
    .trim = FALSE
  )
}

log_info('Writing to {args$output}')
for (run in c('base', sort(unique(mapping$run)))) {
  output_dir <- file.path(args$output, run)
  log_trace('Writing {output_dir}')
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  variables <- list(
    oco2_soundings = args$oco2_soundings,
    tccon_sounding_directory = args$tccon_sounding_directory,
    obspack_directory = args$obspack_directory,
    lauder = args$lauder,
    input_run = file.path(args$input, run),
    meteorology_run = file.path(args$input, 'base'),
    output = output_dir,
    utils_partial = Sys.getenv('UTILS_PARTIAL')
  )
  for (filename in list.files(args$template_directory)) {
    file_text <- read_template_file(filename, variables)
    cat(file_text, file = file.path(output_dir, filename))
  }
}

log_info('Done')
