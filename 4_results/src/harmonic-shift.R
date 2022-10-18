library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(patchwork)
library(lubridate, warn.conflicts = FALSE)
library(tidyr)
library(sf)

source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--climatology-by-region')
parser$add_argument('--regions')
parser$add_argument('--region-sf')
parser$add_argument('--sib4-climatology-assim')
parser$add_argument('--output')
args <- parser$parse_args()

climatology_by_region <- readRDS(args$climatology_by_region) %>%
  filter(inventory != 'ocean')

if (args$regions == 'global-zonal') {
  position_point_dodge <- position_dodge(width = 0.6)
  climatology_by_region <- climatology_by_region %>%
    filter(region == 'Global' | grepl('°', as.character(region)))
  axis_text_theme <- element_text(angle = 15, hjust = 1)
  has_inset <- FALSE
} else if (args$regions == 'land-transcoms') {
  position_point_dodge <- position_dodge(width = 0.8)
  climatology_by_region <- climatology_by_region %>%
    filter(region %in% sprintf('T%02d', 1 : 11))
  axis_text_theme <- element_text(angle = 0, hjust = 0.5)

  has_inset <- TRUE
  region_sf <- readRDS(args$region_sf)
  region_sf_subset <- region_sf %>%
    filter(region_code %in% sprintf('T%02d', 1 : 11))
  region_colours <- region_sf_subset$colour
  names(region_colours) <- region_sf_subset$region_code
  inset_map <- ggplot() +
    geom_sf(
      data = region_sf_subset,
      mapping = aes(fill = region_code),
      colour = '#999999',
      size = 0.1
    ) +
    geom_sf_text(
      data = region_sf_subset,
      aes(label = region_code),
      colour = 'black',
      size = 2,
      nudge_x = region_sf_subset$label_nudge_x,
      nudge_y = region_sf_subset$label_nudge_y
    ) +
    scale_fill_manual(
      values = region_colours
    ) +
    guides(fill = 'none') +
    labs(x = NULL, y = NULL) +
    scale_x_continuous(breaks = c(-179.9, 179.9)) +
    scale_y_continuous(breaks = c(-89.9, 89.9)) +
    coord_sf(
      xlim = c(-180, 180),
      ylim = c(-90, 90),
      default_crs = st_crs('WGS84'),
      label_graticule = '',
      crs = st_crs('+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs')
    ) +
    theme(
      panel.border = element_blank(),
      panel.grid = element_line(colour = 'grey20', size = 0.1)
    )
} else {
  stop('regions argument not supported')
}

with_nc_file(list(fn = args$sib4_climatology_assim), {
  origin <- ymd_hms(ncdf4::ncatt_get(fn, 0, 'origin')$value)
})

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

climatology_by_region_harmonic <- climatology_by_region %>%
  filter(!(variable %in% c('intercept', 'trend'))) %>%
  separate(
    variable,
    c('part', 'harmonic', 'coefficient'),
    sep = '_'
  ) %>%
  filter(harmonic == '1') %>%
  group_by(output, inventory, region, part) %>%
  summarise(
    start = t(value_samples[1, ] + trend_start * value_samples[2, ]),
    end = t(value_samples[1, ] + trend_end * value_samples[2, ])
  ) %>%
  group_by(output, inventory, region) %>%
  summarise(
    amplitude_start = t(sqrt(start[1, ] ^ 2 + start[2, ] ^ 2)),
    amplitude_end = t(sqrt(end[1, ] ^ 2 + end[2, ] ^ 2)),
    phase_start = t(atan(-start[2, ] / start[1, ])),
    phase_end = t(atan(-end[2, ] / end[1, ])),
    groups = '.drop'
  ) %>%
  mutate(
    amplitude_diff = amplitude_end - amplitude_start,
    amplitude_percent = ifelse(
      amplitude_start == 0,
      0,
      100 * (amplitude_end - amplitude_start) / amplitude_start
    ),
    phase_raw_diff = phase_end - phase_start,
    phase_diff = (phase_raw_diff + pi / 2) %% pi - pi / 2,
    phase_diff_days = -ifelse(
      amplitude_start == 0,
      0,
      365.25 * phase_diff / (2 * pi)
    )
  ) %>%
  mutate(
    phase_diff_days_q025 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.025),
    phase_diff_days_q500 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.500),
    phase_diff_days_q975 = matrixStats::rowQuantiles(phase_diff_days, probs = 0.975),
    amplitude_percent_q025 = matrixStats::rowQuantiles(amplitude_percent, probs = 0.025),
    amplitude_percent_q500 = matrixStats::rowQuantiles(amplitude_percent, probs = 0.500),
    amplitude_percent_q975 = matrixStats::rowQuantiles(amplitude_percent, probs = 0.975),
    amplitude_diff_q025 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.025),
    amplitude_diff_q500 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.500),
    amplitude_diff_q975 = matrixStats::rowQuantiles(amplitude_diff, probs = 0.975)
  ) %>%
  mutate(
    region = droplevels(region),
    inventory = factor(c(
      'bio_assim' = 'GPP',
      'bio_resp_tot' = 'Respiration',
      'nee' = 'NEE'
    )[inventory], c(
      'GPP',
      'Respiration',
      'NEE'
    )),
    output = factor(c(
      'Bottom up' = 'Bottom up',
      'Posterior' = 'Posterior\n(95% interval)'
    )[output])
  )

regions <- factor(
  levels(climatology_by_region_harmonic$region),
  levels(climatology_by_region_harmonic$region),
)

plot_changes <- function(...) {
  ggplot(...) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_tile(
      data = data.frame(
        region = regions[seq(2, length(regions), by = 2)]
      ),
      mapping = aes(
        x = region,
        y = 0,
        height = Inf,
        width = 1
      ),
      fill = '#444444',
      alpha = 0.1,
      inherit.aes = FALSE
    ) +
    expand_limits(x = c(0.5, length(regions) + 0.5)) +
    geom_linerange(
      data = climatology_by_region_harmonic %>% filter(startsWith(as.character(output), 'Posterior')),
      position = position_point_dodge,
      size = 0.6
    ) +
    geom_point(
      data = climatology_by_region_harmonic %>% filter(startsWith(as.character(output), 'Bottom')),
      position = position_point_dodge,
      shape = 4,
      size = 2
    ) +
    geom_point(
      data = climatology_by_region_harmonic,
      mapping = aes(shape = output),
      alpha = 0
    ) +
    scale_x_discrete(
      limits = levels(climatology_by_region_harmonic$region),
      expand = c(0, 0, 0, 0)
    ) +
    scale_colour_manual(
      values = c('#018571', '#a6611a', 'black')
    ) +
    scale_shape_manual(
      values = c(4, charToRaw('|')),
      drop = FALSE
    ) +
    labs(
      x = NULL,
      colour = 'Component',
      shape = 'Estimate'
    ) +
    guides(shape = guide_legend(override.aes = list(size = c(3, 6), alpha = 1))) +
    theme(
      axis.title.y = element_text(margin = margin(0, 10, 0, 0, unit = 'points')),
      axis.text.x = axis_text_theme,
      plot.title = element_text(size = 10),
      plot.margin = margin(t = 0.2, b = 0, l = 0, r = 0, unit = 'cm')
    )
}

phase_shift <- plot_changes(
  mapping = aes(
    region,
    y = phase_diff_days_q500,
    ymin = phase_diff_days_q025,
    ymax = phase_diff_days_q975,
    colour = inventory
  )
) +
  geom_segment(
    data = data.frame(
      y = c(
        seq(-100, -20, by = 10),
        seq(-9, -2, by = 1),
        seq(-0.9, -0.1, by = 0.1),
        seq(0.1, 0.9, by = 0.1),
        seq(2, 9, by = 1),
        seq(20, 100, by = 10)
      )
    ) %>%
      filter(
        y >= 10 * floor(min(climatology_by_region_harmonic$phase_diff_days_q025) / 10),
        y <= 10 * ceiling(max(climatology_by_region_harmonic$phase_diff_days_q975) / 10)
      ),
    mapping = aes(x = 0.5, xend = 0.6, y = y, yend = y),
    inherit.aes = FALSE,
    size = 0.1
  ) +
  scale_y_continuous(
    breaks = c(-10, -1, 0, 1, 10),
    trans = scales::pseudo_log_trans(sigma = 0.1, base = 10)
  ) +
  labs(
    y = expression(Delta*tilde(P)['c,1'](italic(S))*' [days]')
  ) +
  ggtitle('Phase shift from January 2015 to December 2020') +
  inset_element(
    grid::textGrob(
      expression('' %<-% 'later   earlier' %->% ''),
      rot = 90,
      gp = grid::gpar(fontsize = 9, col = '#23373b')
    ),
    left = -0.14,
    bottom = 1,
    right = -0.11,
    top = 0,
    align_to = 'panel'
  )

amplitude_change <- plot_changes(
  mapping = aes(
    region,
    y = amplitude_diff_q500,
    ymin = amplitude_diff_q025,
    ymax = amplitude_diff_q975,
    colour = inventory
  )
) +
  labs(
    y = expression(Delta*tilde(A)['c,1'](italic(S))*' [PgC/month]')
  ) +
  ggtitle('Amplitude change from January 2015 to December 2020')

output <- wrap_plots(
  amplitude_change,
  phase_shift,
  ncol = 1,
  guides = 'collect'
)

if (has_inset) {
  output <- output +
    inset_element(
      inset_map,
      left = 0.8,
      bottom = 0,
      right = 1,
      top = 0.7,
      align_to = 'full'
    ) &
    theme(legend.justification = 'top')
}

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$full_width,
  height = 9
)
