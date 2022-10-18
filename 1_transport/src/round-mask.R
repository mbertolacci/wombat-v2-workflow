library(argparse)
library(ncdf4)
library(wombatbasis)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--input')
parser$add_argument('--output')
args <- parser$parse_args()

grid3 <- read_nc_grid_data(
  args$input,
  sprintf('Region%02d', 1 : 22)
)

grid1 <- grid

total_mask <- NULL
for (region_name in names(grid)) {
  value <- round(grid[[region_name]])

  if (is.null(total_mask)) {
    total_mask <- grid[[region_name]]
  } else {
    value <- 0 + (value & !total_mask)
    total_mask <- total_mask | value
  }
  grid[[region_name]] <- value
}

stopifnot(all(range(Reduce(`+`, grid)) == c(0, 1)))

dimensions <- list(
  ncdim_def(
    'lon',
    units = 'degrees_east',
    vals = attr(grid[[1]], 'longitude'),
    longname = 'Longitude'
  ),
  ncdim_def(
    'lat',
    units = 'degrees_north',
    vals = attr(grid[[1]], 'latitude'),
    longname = 'Latitude'
  )
)

output_variables <- lapply(names(grid), function(region_name) {
  ncvar_def(
    region_name,
    units = '1',
    longname = sprintf('Mask for %s', region_name),
    dim = dimensions,
    prec = 'double',
    compression = 6
  )
})

fn <- nc_create(args$output, output_variables)
for (i in seq_along(output_variables)) {
  ncvar_put(fn, output_variables[[i]], grid[[names(grid)[i]]])
}
nc_close(fn)

regions <- 0 : 22
for (i in regions) {
  parts <- sprintf(
    'Region%02d %s Region%02d',
    i,
    ifelse(
      regions >= i,
      '>=',
      '> '
    ),
    regions
  )
  cat(sprintf('Region%02d = Region%02d > 0 && %s;\n', i, i, paste0(parts, collapse = ' && ')))
}
