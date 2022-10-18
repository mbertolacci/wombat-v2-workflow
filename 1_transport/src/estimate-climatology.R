library(argparse)
library(lubridate, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('CLIMATOLOGY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--input')
parser$add_argument('--field-name')
parser$add_argument('--time-resolution')
parser$add_argument('--time-alignment')
parser$add_argument('--origin')
parser$add_argument('--harmonics', type = 'integer')
parser$add_argument('--output')
args <- parser$parse_args()

estimate_climatology(
  args$input,
  args$field_name,
  args$time_resolution,
  args$time_alignment,
  ymd_hms(args$origin),
  args$harmonics,
  args$output
)
