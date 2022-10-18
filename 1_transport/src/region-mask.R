library(argparse)
library(wombatbasis)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--transcom-mask')
parser$add_argument('--output')
args <- parser$parse_args()

log_info('Reading regions from {args$transcom_mask}')
region_grid <- read_nc_grid_data(
  args$transcom_mask,
  sprintf('Region%02d', 1 : 22)
)

# Add the New Zealand land region
nz_region <- region_grid$Region01
nz_region[1 : length(nz_region)] <- 0
nz_region[140, 23] <- 1
nz_region[141, 23] <- 1
nz_region[141, 24] <- 1
nz_region[142, 25] <- 1
nz_region[143, 26] <- 1
nz_region[143, 27] <- 1
nz_region[144, 27] <- 1
for (region_name in names(region_grid)) {
  value <- region_grid[[region_name]]
  value <- 0 + (value & !nz_region)
  attributes(value) <- attributes(nz_region)

  region_grid[[region_name]] <- value
}
region_grid$RegionNZ <- nz_region

# Move the Region10 NZ bits to Region15, an ocean region
region_grid$Region10[142, 24] <- 0
region_grid$Region15[142, 24] <- 1
region_grid$Region10[144, 26] <- 0
region_grid$Region15[144, 26] <- 1
region_grid$Region10[144, 26] <- 0
region_grid$Region15[144, 26] <- 1

# Check every grid cell has only one region
stopifnot(all(range(Reduce(`+`, region_grid)) == c(0, 1)))

log_info('Writing region file to {args$output}')
saveRDS(region_grid, args$output)

log_info('Done')
