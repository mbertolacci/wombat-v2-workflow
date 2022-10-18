source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(patchwork)

geom_world_sf <- function(
  colour = 'black',
  size = 0.1,
  inherit.aes = FALSE,
  ...
) {
  worldmap <- sf::st_wrap_dateline(sf::st_as_sf(
    rnaturalearth::ne_coastline(110, returnclass = 'sf')
  ))

  geom_sf(
    data = worldmap,
    colour = colour,
    size = size,
    inherit.aes = inherit.aes,
    ...
  )
}

match_grid <- function(df, model) {
  latitude_index <- findInterval(
    df$latitude,
    c(model$latitude[1] - model$cell_height[1] / 2, model$latitude + model$cell_height / 2),
    rightmost.closed = TRUE
  )
  longitude_index <- findInterval(
    df$longitude,
    c(model$longitude[1] - model$cell_width[1] / 2, model$longitude + model$cell_width / 2),
    rightmost.closed = TRUE
  )
  longitude_index[longitude_index > length(model$longitude)] <- 1L

  df %>%
    mutate(
      latitude = model$latitude[latitude_index],
      longitude = model$longitude[longitude_index]
    )
}


parser <- ArgumentParser()
parser$add_argument('--observations')
parser$add_argument('--output')
args <- parser$parse_args()

observations <- fst::read_fst(args$observations)

# grid_system <- list(
#   latitude = seq(-85, 85, by = 10),
#   cell_height = 10,
#   longitude = seq(-175, 175, by = 10),
#   cell_width = 10
# )
grid_system <- list(
  latitude = seq(-87.5, 87.5, by = 5),
  cell_height = 5,
  longitude = seq(-177.5, 177.5, by = 5),
  cell_width = 5
)
grid_base <- expand.grid(
  latitude = grid_system$latitude,
  longitude = grid_system$longitude
)

count_df <- lapply(c('oco2', 'obspack'), function(observation_type_i) {
  grid_base %>%
    left_join(
      observations %>%
        filter(
          observation_type == observation_type_i,
          assimilate == 1,
          overall_observation_mode != 'OG'
        ) %>%
        match_grid(grid_system) %>%
        group_by(longitude, latitude) %>%
        summarise(n = n(), .groups = 'drop'),
      by = c('longitude', 'latitude')
    ) %>%
    grid_df_to_sf('n') %>%
    mutate(observation_type = observation_type_i)
}) %>%
  bind_rows() %>%
  mutate(
    observation_type = factor(c(
      'oco2' = 'OCO-2 Land',
      'obspack' = 'In situ/flask'
    )[as.character(observation_type)], c(
      'OCO-2 Land',
      'In situ/flask'
    ))
  )

output1 <- ggplot() +
  geom_sf(
    data = count_df,
    mapping = aes(fill = pmax(10, n)),
    colour = NA
  ) +
  geom_path(
    data = data.frame(
      longitude = c(-180, 180, 180, -180, -180),
      latitude = c(-90, -90, 90, 90, -90)
    ),
    mapping = aes(longitude, latitude),
    colour = 'black',
    size = 0.2
  ) +
  geom_world_sf() +
  scale_fill_binned(
    type = 'viridis',
    trans = 'log10',
    breaks = 10 ^ seq(0, 6, by = 0.5),
    labels = function(x) {
      ifelse(x %% 1 == 0, sprintf('%.0f', x), '')
    },
    na.value = NA
  ) +
  coord_sf(
    crs = sf::st_crs('+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs'),
    default_crs = sf::st_crs('WGS84'),
    label_graticule = '',
    expand = FALSE
  ) +
  labs(x = NULL, y = NULL, fill = '# Observations') +
  facet_wrap(~ observation_type) +
  theme(
    panel.border = element_blank(),
    legend.key.height = unit(0.5, 'cm')
  )

output2 <- observations %>%
  filter(
    assimilate == 1,
    overall_observation_mode != 'OG'
  ) %>%
  mutate(
    month = lubridate::round_date(time, 'month')
  ) %>%
  group_by(observation_type, month) %>%
  summarise(n = n(), .groups = 'drop') %>%
  mutate(
    observation_type = factor(c(
      'oco2' = 'OCO-2 Land',
      'obspack' = 'In situ/flask'
    )[as.character(observation_type)], c(
      'OCO-2 Land',
      'In situ/flask'
    ))
  ) %>%
  ggplot(aes(month, n, linetype = observation_type)) +
    geom_line() +
    labs(x = 'Month', y = '# Observations', linetype = NULL) +
    theme(legend.position = 'right')

output <- wrap_plots(output1, output2, ncol = 1, heights = c(1, 0.6))

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$full_width,
  height = 7
)
