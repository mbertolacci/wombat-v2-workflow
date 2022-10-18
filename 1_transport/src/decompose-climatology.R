library(argparse)
library(lubridate, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('CLIMATOLOGY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--input', nargs = '+')
parser$add_argument('--output', nargs = '+')
parser$add_argument('--climatology')
parser$add_argument('--field-name')
parser$add_argument('--time-resolution')
parser$add_argument('--time-alignment')
args <- parser$parse_args()

decompose_climatology(
  args$input,
  args$output,
  args$climatology,
  args$field_name,
  args$time_resolution,
  args$time_alignment,
  fields = 'residual'
)
