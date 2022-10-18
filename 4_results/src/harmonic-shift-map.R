library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(lubridate, warn.conflicts = FALSE)
library(Matrix)
library(patchwork)
library(tidyr, warn.conflicts = FALSE)
library(sf)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--control-emissions')
parser$add_argument('--inventory-harmonic')
parser$add_argument('--region-sf')
parser$add_argument('--region-grid')
parser$add_argument('--samples')
parser$add_argument('--sib4-climatology-assim')
parser$add_argument('--sib4-climatology-resp-tot')
parser$add_argument('--output')
args <- parser$parse_args()

control_emissions <- fst::read_fst(args$control_emissions)
samples <- readRDS(args$samples)
region_sf <- readRDS(args$region_sf)

match <- stringr::str_match(args$inventory_harmonic, '(\\w+)(-(\\d))?')
args$inventory <- c(
  'gpp' = 'bio_assim',
  'resp' = 'bio_resp_tot',
  'nee' = 'nee'
)[match[2]]
if (is.na(match[4])) {
  args$harmonic <- 1
} else {
  args$harmonic <- as.integer(match[4])
}
args$in_supplement <- !(args$harmonic == 1 && args$inventory == 'nee')

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

X_fit_region_to_2x25 <- with(
  fit_region_grid_df %>%
    mutate(region = factor(region)),
  sparseMatrix(
    i = seq_along(longitude),
    j = as.integer(region),
    x = rep(1, length(longitude)),
    dims = c(length(longitude), nlevels(region))
  )
)[, 2 : 12]

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
    mutate(inventory = 'bio_resp_tot')
)

alpha_modified <- samples$alpha_df %>%
  select(
    inventory,
    region,
    variable = component,
    alpha_mean = value,
    alpha_samples = value_samples
  ) %>%
  mutate(
    variable = factor(as.character(variable), levels(climatology_prior$variable)),
    alpha_samples = alpha_samples[, floor(seq(1, ncol(alpha_samples), length.out = 100))]
  )

climatology_posterior <- climatology_prior %>%
  group_by(variable, inventory) %>%
  group_modify(~ {
    alpha_samples <- alpha_modified %>%
      filter(variable == .y$variable, inventory == .y$inventory) %>%
      pull(alpha_samples)

    if (nrow(alpha_samples) == 0) {
      alpha_samples <- matrix(0, nrow = 11, ncol = ncol(alpha_samples))
    }

    X_i <- Diagonal(x = .x$value) %*% X_fit_region_to_2x25

    .x %>%
      select(-value) %>%
      mutate(
        value_samples = .x$value + as.matrix(X_i %*% alpha_samples),
        value = rowMeans(value_samples)
      )
  })

climatology <- bind_rows(
  climatology_prior %>%
    mutate(
      output = 'Bottom-up',
      value_samples = matrix(
        value,
        nrow = length(value),
        ncol = ncol(climatology_posterior$value_samples)
      )
    ),
  climatology_posterior %>% mutate(output = 'Posterior')
) %>%
  arrange(longitude, latitude) %>%
  {
    x <- .
    bind_rows(
      x,
      x %>%
        group_by(longitude, latitude, output, variable) %>%
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

times <- seq(as.Date('2015-01-01'), as.Date('2020-12-01'), by = '1 month')
X_climatology <- climatology_design_matrix(
  time = times,
  origin = origin,
  time_alignment = 'start',
  time_resolution = 'monthly',
  harmonics = 3
)

trend_start <- min(X_climatology[, 2])
trend_end <- max(X_climatology[, 2])

climatology_start_end <- climatology %>%
  filter(!(variable %in% c('intercept', 'trend'))) %>%
  separate(
    variable,
    c('part', 'harmonic', 'coefficient'),
    sep = '_'
  ) %>%
  mutate(harmonic = as.character(harmonic)) %>%
  select(-value) %>%
  pivot_wider(names_from = 'coefficient', values_from = 'value_samples') %>%
  mutate(
    start = intercept + trend_start * trend,
    end = intercept + trend_end * trend
  ) %>%
  select(-intercept, -trend)

climatology_harmonic_start <- climatology_start_end %>%
  select(-end) %>%
  pivot_wider(names_from = 'part', values_from = 'start') %>%
  mutate(
    amplitude_start = sqrt(cos12 ^ 2 + sin12 ^ 2),
    phase_start = ifelse(cos12 == 0, 0, atan(-sin12 / cos12))
  ) %>%
  select(-cos12, -sin12)
climatology_harmonic_end <- climatology_start_end %>%
  select(-start) %>%
  pivot_wider(names_from = 'part', values_from = 'end') %>%
  mutate(
    amplitude_end = sqrt(cos12 ^ 2 + sin12 ^ 2),
    phase_end = ifelse(cos12 == 0, 0, atan(-sin12 / cos12))
  ) %>%
  select(-cos12, -sin12)

climatology_harmonic <- climatology_harmonic_start
climatology_harmonic$amplitude_end <- climatology_harmonic_end$amplitude_end
climatology_harmonic$phase_end <- climatology_harmonic_end$phase_end

climatology_harmonic_subset <- climatology_harmonic %>%
  filter(
    inventory == args$inventory,
    harmonic == args$harmonic,
    abs(latitude) != 89.5
  ) %>%
  mutate(
    amplitude_diff = amplitude_end - amplitude_start,
    phase_diff = (phase_end - phase_start + pi / 2) %% pi - pi / 2,
    phase_diff_days = -ifelse(
      amplitude_start == 0,
      0,
      365.25 * phase_diff / (2 * pi)
    ),
    phase_diff_days_q250 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.250),
    phase_diff_days_q500 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.500),
    phase_diff_days_q750 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.750),
    phase_diff_days_iqr = phase_diff_days_q750 - phase_diff_days_q250,
    amplitude_diff_q250 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.250),
    amplitude_diff_q500 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.500),
    amplitude_diff_q750 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.750),
    amplitude_diff_iqr = amplitude_diff_q750 - amplitude_diff_q250,
    amplitude_diff_iqr = PER_SECONDS_TO_PER_YEAR * ifelse(amplitude_diff_q500 == 0, NA, amplitude_diff_iqr),
    amplitude_diff_q500 = PER_SECONDS_TO_PER_YEAR * ifelse(amplitude_diff_q500 == 0, NA, amplitude_diff_q500),
    phase_diff_days_iqr = ifelse(phase_diff_days_q500 == 0, NA, phase_diff_days_iqr),
    phase_diff_days_q500 = ifelse(phase_diff_days_q500 == 0, NA, phase_diff_days_q500)
  )

phase_breaks <- c(-5 : 5)
phase_limits <- c(-6.5, 6.5)
phase_iqr_breaks <- phase_breaks[phase_breaks > 0]
phase_iqr_limits <- c(0, phase_limits[2])

if (args$inventory == 'nee') {
  phase_shift_label <- labs(
    fill = expression('' %<-% 'later       '*Delta*P['nee,1'](bold(s))*' [days]       earlier' %->% '')
  )
  phase_shift_iqr_label <- labs(fill = expression('Posterior IQR of '*Delta*P['nee,1'](bold(s))*' [days]'))
} else if (args$inventory == 'bio_assim') {
  phase_shift_label <- labs(
    fill = expression('' %<-% 'later       '*Delta*P['gpp,1'](bold(s))*' [days]       earlier' %->% '')
  )
  phase_shift_iqr_label <- labs(fill = expression('Posterior IQR of '*Delta*P['gpp,1'](bold(s))*' [days]'))
} else if (args$inventory == 'bio_resp_tot') {
  phase_shift_label <- labs(
    fill = expression('' %<-% 'later       '*Delta*P['resp,1'](bold(s))*' [days]       earlier' %->% '')
  )
  phase_shift_iqr_label <- labs(fill = expression('Posterior IQR of '*Delta*P['resp,1'](bold(s))*' [days]'))
} else {
  stop('unsupported inventory')
}

phase_bottom_up <- climatology_harmonic_subset %>%
  filter(output == 'Bottom-up') %>%
  mutate(
    phase_diff_days_q500 = discretise_by_breaks(phase_diff_days_q500, phase_breaks, phase_limits)
  ) %>%
  grid_df_to_sf('phase_diff_days_q500') %>%
    plot_map(
      phase_diff_days_q500,
      phase_breaks,
      phase_limits,
      drop_second_labels = FALSE,
      show_excess = TRUE,
      label_precision = '0'
    ) +
    phase_shift_label +
    ggtitle('Bottom-up')

phase_posterior <- climatology_harmonic_subset %>%
  filter(output == 'Posterior') %>%
  mutate(
    phase_diff_days_q500 = discretise_by_breaks(phase_diff_days_q500, phase_breaks, phase_limits)
  ) %>%
  grid_df_to_sf('phase_diff_days_q500') %>%
    plot_map(
      phase_diff_days_q500,
      phase_breaks,
      phase_limits,
      drop_second_labels = FALSE,
      show_excess = TRUE,
      label_precision = '0'
    ) +
    phase_shift_label +
    ggtitle('Posterior median')

phase_iqr <- climatology_harmonic_subset %>%
  filter(output == 'Posterior') %>%
  mutate(
    phase_diff_days_iqr = discretise_by_breaks(phase_diff_days_iqr, phase_iqr_breaks, phase_iqr_limits)
  ) %>%
  grid_df_to_sf('phase_diff_days_iqr') %>%
    plot_map(
      phase_diff_days_iqr,
      phase_iqr_breaks,
      phase_iqr_limits,
      symmetric = FALSE,
      drop_second_labels = FALSE,
      show_excess = TRUE,
      label_precision = '0'
    ) +
    phase_shift_iqr_label +
    ggtitle('Posterior IQR')

amplitude_breaks <- seq(-0.25, 0.25, by = 0.05)
amplitude_limits <- c(-0.25, 0.25) + c(-1, 1) * 0.05 * 1.5
amplitude_iqr_breaks <- seq(0.01, 0.05, by = 0.01)
amplitude_iqr_limits <- c(0, 0.065)

if (args$inventory == 'nee') {
  amplitude_label <- labs(fill = expression(Delta*A['nee,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year]'))
  amplitude_iqr_label <- labs(fill = expression('Posterior IQR of '*Delta*A['nee,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year]'))
} else if (args$inventory == 'bio_assim') {
  amplitude_label <- labs(fill = expression(Delta*A['gpp,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year]'))
  amplitude_iqr_label <- labs(fill = expression('Posterior IQR of '*Delta*A['gpp,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year]'))
} else if (args$inventory == 'bio_resp_tot') {
  amplitude_label <- labs(fill = expression(Delta*A['resp,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year]'))
  amplitude_iqr_label <- labs(fill = expression('Posterior IQR of '*Delta*A['resp,1'](bold(s))*' [kgCO'[2]*'/'*m^2*'/year]'))
} else {
  stop('unsupported inventory')
}

amplitude_bottom_up <- climatology_harmonic_subset %>%
  filter(output == 'Bottom-up') %>%
  mutate(
    amplitude_diff_q500 = discretise_by_breaks(amplitude_diff_q500, amplitude_breaks, amplitude_limits)
  ) %>%
  grid_df_to_sf('amplitude_diff_q500') %>%
    plot_map(
      amplitude_diff_q500,
      amplitude_breaks,
      amplitude_limits,
      drop_second_labels = TRUE,
      show_excess = TRUE,
      label_precision = '2'
    ) +
    amplitude_label +
    ggtitle('Bottom-up')

amplitude_posterior <- climatology_harmonic_subset %>%
  filter(output == 'Posterior') %>%
  mutate(
    amplitude_diff_q500 = discretise_by_breaks(amplitude_diff_q500, amplitude_breaks, amplitude_limits)
  ) %>%
  grid_df_to_sf('amplitude_diff_q500') %>%
    plot_map(
      amplitude_diff_q500,
      amplitude_breaks,
      amplitude_limits,
      drop_second_labels = TRUE,
      show_excess = TRUE,
      label_precision = '2'
    ) +
    amplitude_label +
    ggtitle('Posterior median')

amplitude_iqr <- climatology_harmonic_subset %>%
  filter(output == 'Posterior') %>%
  mutate(
    amplitude_diff_iqr = discretise_by_breaks(amplitude_diff_iqr, amplitude_iqr_breaks, amplitude_iqr_limits)
  ) %>%
  grid_df_to_sf('amplitude_diff_iqr') %>%
    plot_map(
      amplitude_diff_iqr,
      amplitude_iqr_breaks,
      amplitude_iqr_limits,
      show_excess = TRUE,
      symmetric = FALSE,
      label_precision = '2'
    ) +
    amplitude_iqr_label +
    ggtitle('Posterior IQR')

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

top_theme <- base_theme +
  theme(
    legend.position = 'none',
    plot.margin = margin(t = 0.3, b = 0.2, l = 0.1, r = 0, unit = 'cm')
  )

iqr_theme <- base_theme +
  theme(
    plot.margin = margin(t = 0.4, b = 0.2, l = 0.1, r = 0, unit = 'cm')
  )

phase_column <- wrap_plots(
  phase_bottom_up + top_theme,
  phase_posterior + base_theme,
  phase_iqr + iqr_theme,
  ncol = 1
)

amplitude_column <- wrap_plots(
  amplitude_bottom_up + top_theme,
  amplitude_posterior + base_theme,
  amplitude_iqr + iqr_theme,
  ncol = 1
)

output <- wrap_plots(
  wrap_plots(
    wrap_elements(
      panel = grid::textGrob(sprintf(
        'Amplitude change of %s from\nJanuary 2015 to December 2020',
        c('nee' = 'NEE', 'bio_assim' = 'GPP', 'bio_resp_tot' = 'respiration')[args$inventory]
      ))
    ),
    amplitude_column,
    heights = c(0.07, 1)
  ),
  wrap_plots(
    wrap_elements(
      panel = grid::textGrob(sprintf(
        'Phase shift in %s from\nJanuary 2015 to December 2020',
        c('nee' = 'NEE', 'bio_assim' = 'GPP', 'bio_resp_tot' = 'respiration')[args$inventory]
      ))
    ),
    phase_column,
    heights = c(0.07, 1)
  ),
  nrow = 1
)

ggsave_base(
  args$output,
  output,
  width = if (args$in_supplement) {
    DISPLAY_SETTINGS$supplement_full_width
  } else {
    DISPLAY_SETTINGS$full_width
  },
  height = if (args$in_supplement) 19.5 else 17.8
)
