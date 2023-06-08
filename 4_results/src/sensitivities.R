source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))
library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(ncdf4)

parser <- ArgumentParser()
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations')
parser$add_argument('--xco2-daily-base')
parser$add_argument('--xco2-daily-201601r02')
parser$add_argument('--xco2-daily-201601r06')
parser$add_argument('--output')
args <- parser$parse_args()

args <- list(
  region_sf = '4_results/intermediates/region-sf.rds',
  control_emissions = 'data/wombat-v1-control-emissions.fst',
  perturbations = 'data/wombat-v1-perturbations.fst',
  xco2_daily_base = 'data/wombat-v1-xco2-daily-base.nc',
  xco2_daily_201601r02 = 'data/wombat-v1-xco2-daily-201601R02.nc',
  xco2_daily_201601r06 = 'data/wombat-v1-xco2-daily-201601R06.nc',
  output = '4_results/figures/sensitivities.pdf'
)

region_sf <- readRDS(args$region_sf)

GEOS_CHEM_GRID <- list(
  longitude = list(
    centres = seq(-180, 177.51, 2.5),
    widths = rep(2.5, 144)
  ),
  latitude = list(
    centres = c(-89.5, seq(-88, 88.1, 2), 89.5),
    widths = c(1, rep(2, 89), 1)
  )
)

clamp <- function(x, a, b) pmin(pmax(x, a), b)

low_colour <- '#35978f'
mid_colour <- '#ffffff'
high_colour <- '#bf812d'

flux_breaks <- seq(-1.5, 1.5, by = 0.25)
flux_limits <- c(-1.5, 1.5)

xco2_breaks <- seq(-0.3, 0.3, by = 0.05)
xco2_limits <- c(-0.3, 0.3)

read_xco2_daily <- function(filename) {
  with_nc_file(list(fn = filename), {
    times <- ymd_hms('2000-01-01 00:00:00') + minutes(ncvar_get(fn, 'time'))
    xco2 <- ncvar_get(fn, 'xco2')
  })

  expand.grid(
    longitude_index = seq_along(GEOS_CHEM_GRID$longitude$centres),
    latitude_index = seq_along(GEOS_CHEM_GRID$latitude$centres),
    time = times
  ) %>%
    mutate(
      date = as.Date(time),
      longitude = GEOS_CHEM_GRID$longitude$centres[longitude_index],
      cell_width = GEOS_CHEM_GRID$longitude$widths[longitude_index],
      latitude = GEOS_CHEM_GRID$latitude$centres[latitude_index],
      cell_height = GEOS_CHEM_GRID$latitude$widths[latitude_index],
      xco2 = as.vector(xco2)
    ) %>%
    select(-time, -longitude_index, -latitude_index)
}

xco2_base_df <- read_xco2_daily(args$xco2_daily_base)
xco2_201601r02_df <- read_xco2_daily(args$xco2_daily_201601r02) %>%
  mutate(region = 2)
xco2_201601r06_df <- read_xco2_daily(args$xco2_daily_201601r06) %>%
  mutate(region = 6)

xco2_sensitivity <- bind_rows(
  xco2_201601r02_df,
  xco2_201601r06_df
) %>%
  left_join(
    xco2_base_df %>% select(date, longitude, latitude, xco2_base = xco2),
    by = c('date', 'longitude', 'latitude')
  ) %>%
  mutate(
    xco2_sensitivity = xco2 - xco2_base
  )

control_emissions <- fst::read_fst(args$control_emissions)
perturbations <- fst::read_fst(args$perturbations) %>%
  left_join(
    control_emissions %>%
      select(
        model_id,
        month_start,
        longitude,
        cell_width,
        latitude,
        cell_height
      ),
    by = 'model_id'
  ) %>%
  mutate(flux_density = flux_density * 31536000)

perturbations_sf <- perturbations %>%
  filter(region %in% c(2, 6), month_start == '2016-01-01') %>%
  mutate(
    date = '2016-01',
    region = sprintf('T%02d', region)
  ) %>%
  group_by(date, region) %>%
  group_map(~ {
    .x %>%
      filter(
        abs(latitude) != 89.5
      ) %>%
      arrange(longitude, latitude) %>%
      mutate(
        value = discretise_by_breaks(
          ifelse(abs(flux_density) <= 0.01, NA, flux_density),
          flux_breaks,
          flux_limits
        )
      ) %>%
      grid_df_to_sf('value') %>%
      mutate(
        date = .y$date,
        region = .y$region
      )
  }) %>%
  bind_rows() %>%
  mutate(region_name = factor(region))

sensitivities_sf <- xco2_sensitivity %>%
  filter(
    region %in% c(2, 6),
    format(date) %in% c(
      '2016-01-01',
      '2016-01-15',
      '2016-02-15'
    )
  ) %>%
  mutate(
    region = sprintf('T%02d', region)
  ) %>%
  group_by(date, region) %>%
  group_map(~ {
    .x %>%
      filter(
        abs(latitude) != 89.5
      ) %>%
      arrange(longitude, latitude) %>%
      mutate(
        value = discretise_by_breaks(
          ifelse(
            abs(xco2_sensitivity) < 0.01,
            NA,
            ifelse(
              # HACK(mgnb): remove weird numerical artefacts
              format(.y$date) == '2016-01-01'
              & (
                (
                  .y$region == 'T02'
                  & (latitude < 15 | longitude > 0)
                )
                | (
                  .y$region == 'T06'
                  & (latitude > 10 | longitude < 0)
                )
              ),
              NA,
              xco2_sensitivity
            )
          ),
          xco2_breaks,
          xco2_limits
        )
      ) %>%
      grid_df_to_sf('value') %>%
      mutate(
        date = .y$date,
        region = .y$region
      )
  }) %>%
  bind_rows() %>%
  mutate(region_name = factor(region))

flux_plot <- plot_map(
  perturbations_sf,
  value,
  flux_breaks,
  flux_limits,
  show_excess = TRUE,
  label_precision = 2,
  bar_width = 16,
  drop_second_labels = TRUE
) +
  facet_grid(date ~ region_name) +
  labs(fill = expression('Flux [kg/'*m^2*'/year]')) +
  ggtitle(expression('Pulse in flux field, '*X(bold(s), t)))

sensitivities_plot <- plot_map(
  sensitivities_sf,
  value,
  xco2_breaks,
  xco2_limits,
  show_excess = TRUE,
  label_precision = 2,
  bar_width = 16,
  drop_second_labels = TRUE
) +
  facet_grid(date ~ region_name) +
  labs(fill = expression('XCO'[2]*' response [ppm]')) +
  ggtitle(expression('Response in mole-fraction field, '*Y(bold(s), h, t)))

base_theme <- theme(
  legend.position = 'bottom',
  legend.margin = margin(t = -0.2, l = 0, b = -0.2, r = 0, unit = 'cm'),
  legend.title = element_text(size = 9),
  plot.title = element_text(
    hjust = 0.5,
    vjust = 1,
    size = 11 ,
    margin = margin(t = 0, r = 0, b = 0, l = 0, unit = 'cm')
  ),
  plot.margin = margin(t = 0, b = 0.3, l = 0.1, r = 0, unit = 'cm')
)

output <- wrap_plots(
  flux_plot + base_theme,
  sensitivities_plot +
    base_theme +
    theme(
      plot.margin = margin(t = 0.5, b = 0.3, l = 0.1, r = 0, unit = 'cm')
    ),
  ncol = 1,
  heights = c(1, 3.13)
)

ggsave(
  args$output,
  plot = output,
  width = DISPLAY_SETTINGS$full_width,
  height = 20,
  units = 'cm'
)
