source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(sf)

parser <- ArgumentParser()
parser$add_argument('--region-sf')
parser$add_argument('--output')
args <- parser$parse_args()

region_sf <- readRDS(args$region_sf)

region_colours <- region_sf$colour
names(region_colours) <- region_sf$region_code

output <- ggplot() +
  geom_sf(
    data = region_sf,
    mapping = aes(fill = region_code),
    colour = '#555555',
    size = 0.1
  ) +
  geom_sf_text(
    data = region_sf,
    mapping = aes(label = region_code),
    colour = 'black',
    nudge_x = region_sf$label_nudge_x,
    nudge_y = region_sf$label_nudge_y
  ) +
  scale_fill_manual(values = region_colours) +
  guides(fill = 'none') +
  labs(x = NULL, y = NULL) +
  coord_sf(
    default_crs = sf::st_crs('WGS84'),
    crs = st_crs('+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs')
  ) +
  theme(
    panel.border = element_blank(),
    panel.grid = element_line(colour = 'grey20', size = 0.1) #,
  )

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 8.5
)
