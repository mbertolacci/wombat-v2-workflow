library(ggplot2)

theme_set(theme_bw())
theme_replace(
  legend.background = element_blank(),
  legend.key = element_blank(),
  panel.background = element_blank(),
  strip.background = element_blank(),
  plot.background = element_blank(),
  # panel.border = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  axis.text = element_text(colour = '#23373b'),
  axis.title = element_text(colour = '#23373b'),
  plot.margin = margin(t = 1, r = 0, b = 0, l = 1, unit = 'mm')
)

DISPLAY_SETTINGS <- list(
  full_width = 14.32,
  supplement_full_width = 16.5,
  full_height = 20,
  png_plot_dpi = 320
)

get_legend <- function(plot_object){
  tmp <- ggplot_gtable(ggplot_build(plot_object))
  legend_index <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
  tmp$grobs[[legend_index]]
}

ggsave_base <- function(filename, plot, bg = 'transparent', dpi = DISPLAY_SETTINGS$png_plot_dpi, ...) {
  ggsave(
    filename,
    plot,
    units = 'cm',
    dpi = dpi,
    bg = bg,
    ...
  )
}

ggsave_fullwidth <- function(...) {
  ggsave_base(..., width = DISPLAY_SETTINGS$full_width)
}

scale_fill_fermenter_n <- function(
  ...,
  type = 'seq',
  palette = 1,
  direction = -1,
  n_colours = 10,
  na.value = 'grey50',
  guide = 'coloursteps',
  aesthetics = 'fill'
) {
  palette <- scales::gradient_n_pal(
    scales::brewer_pal(type, palette, direction)(n_colours),
    NULL,
    'Lab'
  )
  binned_scale(
    aesthetics,
    'steps2_uneven',
    function(x) {
      palette(seq(0, 1, length.out = length(x)))
    },
    na.value = na.value,
    guide = guide,
    ...
  )
}

plot_map <- function(
  df,
  variable,
  breaks,
  limits,
  show_excess = TRUE,
  label_precision = '0',
  drop_second_labels = FALSE,
  symmetric = TRUE,
  bar_width = 13
) {
  base_labels <- sprintf(paste0('%.', label_precision, 'f'), breaks)
  labels <- if (show_excess) {
    ifelse(
      abs(breaks) == max(abs(breaks)),
      sprintf(
        paste0('%s%.', label_precision, 'f'),
        ifelse(breaks < 0, '<-', '>'),
        max(breaks)
      ),
      base_labels
    )
  } else {
    base_labels
  }
  if (drop_second_labels) {
    labels[seq(2, length(breaks), by = 2)] <- ''
  }
  ggplot() +
    geom_sf(
      data = df,
      mapping = aes(fill = {{ variable }}),
      color = NA
    ) +
    geom_sf(data = region_sf, fill = NA, colour = '#888888', size = 0.1) +
    geom_segment(
      data = data.frame(y = c(-23, 23, 50)),
      mapping = aes(x = -180, y = y, xend = 180, yend = y),
      colour = 'black',
      linetype = 'dashed',
      size = 0.4
    ) +
    geom_text(
      data = data.frame(
        x = c(-170, -180, -180),
        y = c(-23, 23, 50),
        label = c('23°S', '23°N', '50°N')
      ),
      mapping = aes(x = x, y = y, label = label),
      nudge_y = 10
    ) +
    scale_fill_fermenter_n(
      breaks = breaks,
      palette = if (symmetric) 'RdYlBu' else 'PuBu',
      direction = if (symmetric) -1 else 1,
      n_colours = if (symmetric) 10 else 9,
      limits = limits,
      labels = labels,
      guide = guide_coloursteps(
        title.position = 'top',
        title.hjust = 0.5,
        axis = FALSE,
        label.theme = element_text(size = 8),
        frame.colour = '#999999',
        barwidth = bar_width,
        even.steps = FALSE
      ),
      na.value = '#cccccc'
    ) +
    coord_sf(
      default_crs = sf::st_crs('WGS84'),
      crs = sf::st_crs('+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs')
    ) +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank()
    ) +
    labs(x = NULL, y = NULL)
}
