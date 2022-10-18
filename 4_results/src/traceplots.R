source(Sys.getenv('UTILS_PARTIAL'))
source(Sys.getenv('DISPLAY_PARTIAL'))

library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(patchwork)
library(ggrepel)

parser <- ArgumentParser()
parser$add_argument('--samples')
parser$add_argument('--output')
args <- parser$parse_args()

samples <- readRDS(args$samples)

plot_traces <- function(x, names) {
  if (!is.matrix(x)) x <- t(t(x))
  iterations <- (N_MCMC_WARM_UP + 1) : N_MCMC_SAMPLES
  df <- data.frame(
    iteration = rep(iterations, ncol(x)),
    value = as.vector(x[iterations, ])
  )

  if (!missing(names)) {
    df$name <- rep(names, each = N_MCMC_SAMPLES - N_MCMC_WARM_UP)

    output <- ggplot(mapping = aes(iteration, value, colour = name))
  } else {
    output <- ggplot(mapping = aes(iteration, value))
  }

  output <- output +
    geom_line(data = df)

  if (!missing(names)) {
    label_iteration <- round((N_MCMC_WARM_UP + 1 + N_MCMC_SAMPLES) / 2)
    df <- data.frame(
      name = names,
      iteration = label_iteration,
      value = colMeans(x[iterations, ])
    )
    output <- output + if (length(names) > 2) {
      geom_label_repel(
        data = df,
        size = 3,
        mapping = aes(label = name),
        direction = 'x',
        seed = 3
      )
    } else {
      geom_label(
        data = df,
        size = 3,
        mapping = aes(label = name)
      )
    }
  }

  output +
    guides(colour = 'none') +
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    labs(x = 'Iteration', y = 'Value', colour = NULL) +
    theme(
      plot.margin = margin(t = 1, r = 2, b = 0, l = 0, unit = 'mm')
    )

}

to_bio_labels <- function(x) {
  c('bio_assim' = 'gpp', 'bio_resp_tot' = 'resp')[colnames(x)]
}

scale_bio_inventory <- scale_colour_manual(
  values = c('#018571', '#a6611a')
)

output <- wrap_plots(
  plot_traces(
    samples$alpha[, c('bio_assim.cos12_1_intercept.Region01', 'bio_resp_tot.cos12_1_intercept.Region01')],
    c('gpp', 'resp')
  ) +
    ggtitle(expression(alpha['c,2,1'])) +
    scale_bio_inventory,
  plot_traces(
    samples$alpha[, c('bio_assim.residual.Region05.2015-01', 'bio_resp_tot.residual.Region05.2015-01')],
    c('gpp', 'resp')
  ) +
    ggtitle(expression(alpha['c,6,5,5'])) +
    scale_bio_inventory,
  plot_traces(samples$w_bio_clim, to_bio_labels(samples$w_bio_clim)) +
    ggtitle(expression(tau[c]^beta)) +
    scale_bio_inventory,
  plot_traces(samples$rho_bio_clim) +
    ggtitle(expression(rho['gpp,resp']^beta)),
  plot_traces(samples$w_bio_resid, to_bio_labels(samples$w_bio_resid)) +
    ggtitle(expression(tau[c]^epsilon)) +
    scale_bio_inventory,
  plot_traces(samples$rho_bio_resid) +
    ggtitle(expression(rho['gpp,resp']^epsilon)),
  plot_traces(samples$kappa_bio_resid) +
    ggtitle(expression(kappa['bio']^epsilon)),
  plot_traces(
    samples$gamma,
    c(
      '1_LNLG' = 'OCO-2 land',
      'shipboard' = 'Shipboard',
      'tower' = 'Tower',
      'surface' = 'Surface',
      'aircraft' = 'Aircraft'
    )[colnames(samples$gamma)]
  ) +
    scale_colour_brewer(palette = 'Dark2') +
    ggtitle(expression(gamma[g]^Z)),
  ncol = 2
)

ggsave_base(
  args$output,
  output,
  width = DISPLAY_SETTINGS$supplement_full_width,
  height = 15
)
