library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(Matrix)

source(Sys.getenv('UTILS_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--control-emissions')
parser$add_argument('--region-grid')
parser$add_argument('--samples')
parser$add_argument('--sib4-climatology-assim')
parser$add_argument('--sib4-climatology-resp-tot')
parser$add_argument('--landschutzer-climatology')
parser$add_argument('--output')
args <- parser$parse_args()

control_emissions <- fst::read_fst(args$control_emissions)
samples <- readRDS(args$samples)

region_grid <- readRDS(args$region_grid)
fit_region_grid_df <- NULL
for (latitude_index in seq_along(attr(region_grid$Region01, 'latitude'))) {
  for (longitude_index in seq_along(attr(region_grid$Region01, 'longitude'))) {
    found_any <- FALSE
    for (region_name in names(region_grid)) {
      if (region_grid[[region_name]][longitude_index, latitude_index] > 0) {
        fit_region_grid_df <- rbind(fit_region_grid_df, data.frame(
          longitude = attr(region_grid$Region01, 'longitude')[longitude_index],
          latitude = attr(region_grid$Region01, 'latitude')[latitude_index],
          region = region_name
        ))
        found_any <- TRUE
      }
    }
    if (!found_any) {
      fit_region_grid_df <- rbind(fit_region_grid_df, data.frame(
        longitude = attr(region_grid$Region01, 'longitude')[longitude_index],
        latitude = attr(region_grid$Region01, 'latitude')[latitude_index],
        region = 'Region00'
      ))
    }
  }
}

fit_region_grid_df <- fit_region_grid_df %>%
  mutate(region = factor(region)) %>%
  left_join(
    control_emissions %>%
      distinct(longitude, latitude, area),
    by = c('longitude', 'latitude')
  )

fit_region_df <- fit_region_grid_df %>%
  mutate(
    region = factor(
      ifelse(
        region == 'RegionNZ',
        'NZ',
        stringr::str_replace(as.character(region), 'Region', 'T')
      ),
      levels = c(
        sprintf('T%02d', 0 : 22),
        'NZ'
      )
    )
  ) %>%
  distinct(region) %>%
  arrange(region)
X_2x25_to_fit_region <- with(fit_region_grid_df, sparseMatrix(
  i = as.integer(region),
  j = seq_along(longitude),
  x = area,
  dims = c(nlevels(region), length(longitude))
))

zonal_band_grid_df <- control_emissions %>%
  mutate(
    latitude_bottom = latitude - cell_height / 2
  ) %>%
  distinct(longitude, latitude, latitude_bottom, area) %>%
  arrange(latitude, longitude) %>%
  mutate(
    index = seq_len(n()),
    region = case_when(
      latitude_bottom < -23 ~ '90°S-23°S',
      # NOTE: this is not quite right because there is a grid cell that spans
      # zero latitude; this is fixed below
      latitude_bottom < 0 ~ '23°S-0°',
      latitude_bottom < 23 ~ '0°-23°N',
      latitude_bottom < 55 ~ '23°N-55°N',
      TRUE ~ '55°N-90°N'
    )
  )

zonal_regions <- c(
  '55°N-90°N',
  '23°N-55°N',
  '0°-23°N',
  '23°S-0°',
  '90°S-23°S'
)

zonal_band_grid_with_tropical_split_df <- bind_rows(
  zonal_band_grid_df %>%
    filter(latitude != 0),
  zonal_band_grid_df %>%
    filter(latitude == 0) %>%
    mutate(area = area / 2),
  zonal_band_grid_df %>%
    filter(latitude == 0) %>%
    mutate(
      region = '23°S-0°',
      area = area / 2
    )
) %>%
  mutate(
    region = factor(region, levels = zonal_regions)
  )

zonal_band_df <- zonal_band_grid_with_tropical_split_df %>%
  distinct(region) %>%
  arrange(region)
X_2x25_to_zonal_band <- with(zonal_band_grid_with_tropical_split_df, sparseMatrix(
  i = as.integer(region),
  j = index,
  x = area,
  dims = c(nlevels(region), max(index))
))

global_df <- data.frame(region = 'Global')
X_2x25_to_global <- with(zonal_band_grid_df, sparseMatrix(
  i = rep(1, length(longitude)),
  j = seq_along(longitude),
  x = area,
  dims = c(1, length(longitude))
))

region_df <- bind_rows(global_df, zonal_band_df, fit_region_df) %>%
  mutate(
    region = factor(region, c(
      'Global',
      zonal_regions,
      sprintf('T%02d', 0 : 22),
      'NZ'
    ))
  )
X_2x25_to_regions <- KG_M2_S_TO_PGC_MONTH * rbind(
  X_2x25_to_global,
  X_2x25_to_zonal_band,
  X_2x25_to_fit_region
)

X_fit_region_to_2x25 <- with(fit_region_grid_df, sparseMatrix(
  i = seq_along(longitude),
  j = as.integer(region),
  x = rep(1, length(longitude)),
  dims = c(length(longitude), nlevels(region))
))[, 2 : 12]
X_fit_region_to_regions <- X_2x25_to_regions %*% X_fit_region_to_2x25

origin <- ymd_hms(attr(read_gridded_data(
  args$sib4_climatology_assim,
  'coefficient',
  include_time = FALSE,
  include_variable = TRUE
), 'origin'))

climatology_prior <- bind_rows(
  read_climatology(args$sib4_climatology_assim) %>%
    mutate(
      value = -value,
      inventory = 'bio_assim'
    ),
  read_climatology(args$sib4_climatology_resp_tot) %>%
    mutate(inventory = 'bio_resp_tot'),
  read_climatology(args$landschutzer_climatology) %>%
    mutate(inventory = 'ocean')
)

climatology_prior_region <- climatology_prior %>%
  group_by(variable, inventory) %>%
  group_modify(~ {
    region_df %>%
      mutate(
        value = as.vector(X_2x25_to_regions %*% .x$value)
      )
  })

alpha_modified <- samples$alpha_df %>%
  select(
    inventory,
    region,
    variable = component,
    alpha_mean = value,
    alpha_samples = value_samples
  ) %>%
  mutate(
    variable = factor(as.character(variable), levels(climatology_prior$variable))
  )

climatology_posterior_region <- climatology_prior %>%
  group_by(variable, inventory) %>%
  group_modify(~ {
    alpha_samples <- alpha_modified %>%
      filter(variable == .y$variable, inventory == .y$inventory) %>%
      pull(alpha_samples)

    if (nrow(alpha_samples) == 0) {
      alpha_samples <- matrix(0, nrow = 11, ncol = ncol(alpha_samples))
    }

    value_prior <- as.vector(X_2x25_to_regions %*% .x$value)
    X_i <- X_2x25_to_regions %*% Diagonal(x = .x$value) %*% X_fit_region_to_2x25

    region_df %>%
      mutate(
        value_samples = value_prior + as.matrix(X_i %*% alpha_samples),
        value = rowMeans(value_samples)
      )
  })

climatology_region <- bind_rows(
  climatology_prior_region %>%
    mutate(
      output = 'Bottom up',
      value_samples = matrix(
        value,
        nrow = length(value),
        ncol = ncol(climatology_posterior_region$value_samples)
      )
    ),
  climatology_posterior_region %>% mutate(output = 'Posterior')
) %>%
  {
    x <- .
    bind_rows(
      x,
      x %>%
        filter(inventory %in% c('bio_assim', 'bio_resp_tot')) %>%
        group_by(output, region, variable) %>%
        summarise(
          value = sum(value),
          value_samples = t(colSums(value_samples)),
          .groups = 'drop'
        ) %>%
        mutate(
          inventory = 'nee'
        )
    )
  }

saveRDS(climatology_region, args$output)
