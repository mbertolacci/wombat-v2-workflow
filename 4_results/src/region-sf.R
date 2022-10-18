source(Sys.getenv('UTILS_PARTIAL'))

library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(sf)

parser <- ArgumentParser()
parser$add_argument('--region-grid')
parser$add_argument('--output')
args <- parser$parse_args()

region_grid <- readRDS(args$region_grid)

grid_latitude <- attr(region_grid$Region01, 'latitude')
grid_longitude <- attr(region_grid$Region01, 'longitude')

subset_grid <- function(x) t(x[
  seq_along(grid_longitude),
  (length(grid_latitude) - 1) : 2
])

region_values <- subset_grid(region_grid$Region01)
for (i in 2 : length(region_grid)) {
  x <- subset_grid(region_grid[[i]])
  region_values[x == 1] <- i
}
region_values[region_values == 0] <- NA

# First grid cell is centred on -180, which makes it hard to construct
# compliant polygons; avoid this by splitting horizontal cells in two
region_values_wide <- matrix(
  0,
  nrow = nrow(region_values),
  ncol = 2 * ncol(region_values)
)
for (i in seq_len(nrow(region_values))) {
  # First original cell sits on the boundary, so it's value goes to first and
  # last new cell
  region_values_wide[i, 1] <- region_values[i, 1]
  region_values_wide[i, ncol(region_values_wide)] <- region_values[i, 1]
  # Remaining cells simply repeat
  region_values_wide[
    i,
    2 : (ncol(region_values_wide) - 1)
  ] <- rep(region_values[i, 2 : ncol(region_values)], each = 2)
}

region_raster <- raster::raster(
  region_values_wide,
  xmn = -180,
  xmx = 180,
  ymn = grid_latitude[2] - 1,
  ymx = grid_latitude[length(grid_latitude) - 1] + 1
)
names(region_raster) <- 'region'

region_poly <- region_raster %>%
  raster::rasterToPolygons(dissolve = TRUE)

region_sf <- region_poly %>%
  st_as_sf() %>%
  mutate(
    region_code = ifelse(
      region != 23,
      sprintf('T%02d', region),
      'NZ'
    )
  )
st_crs(region_sf) <- 'WGS84'

ocean_region_colours <- region_sf %>%
  filter(region >= 12, region != 23) %>%
  spdep::poly2nb() %>%
  igraph::graph_from_adj_list() %>%
  igraph::greedy_vertex_coloring()

land_region_nb <- region_sf %>%
  filter(region <= 11 | region == 23) %>%
  spdep::poly2nb()
land_region_nb[[2]] <- integer(0)
land_region_nb[[5]] <- integer(0)

land_region_colours <- land_region_nb %>%
  igraph::graph_from_adj_list() %>%
  igraph::greedy_vertex_coloring()

region_colours <- rep('black', 23)
region_colours[with(region_sf,
  region >= 12 & region != 23
)] <- c(
  '#bbbbff',
  '#9999ff',
  '#7777ff'
)[ocean_region_colours]
region_colours[with(region_sf,
  region <= 11 | region == 23
)] <- c(
  '#bbffbb',
  '#77ff77',
  '#99ff99'
)[land_region_colours]
names(region_colours) <- region_sf$region_code
# region_colours['T10'] <- '#bbffbb'

region_sf$colour <- region_colours

region_sf$label_nudge_x <- 0
region_sf$label_nudge_y <- 0

region_sf$label_nudge_y[region_sf$region_code == 'T08'] <- 5
region_sf$label_nudge_x[region_sf$region_code == 'T09'] <- 15
region_sf$label_nudge_y[region_sf$region_code == 'T09'] <- -10

saveRDS(region_sf, args$output)
