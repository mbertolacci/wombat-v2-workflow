library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Rcpp)

source(Sys.getenv('UTILS_PARTIAL'))

sourceCpp('partials/utils.cpp')

parser <- ArgumentParser()
parser$add_argument('--basis-vectors')
parser$add_argument('--overall-observation-mode', nargs = '+')
parser$add_argument('--control')
parser$add_argument('--observations')
parser$add_argument('--sensitivities-base', nargs = '+')
parser$add_argument('--output')
args <- parser$parse_args()

log_debug('Loading basis vectors from {args$basis_vectors}')
basis_vectors <- fst::read_fst(args$basis_vectors)

log_debug('Loading control from {args$control}')
control <- fst::read_fst(args$control) %>%
  mutate(observation_id = droplevels(observation_id))

log_debug('Loading observations from {args$observations}')
observations <- fst::read_fst(args$observations) %>%
  filter(
    assimilate %in% c(1, 3),
    overall_observation_mode %in% args$overall_observation_mode,
    observation_id %in% control$observation_id
  ) %>%
  arrange(observation_group, time) %>%
  mutate(
    observation_id = factor(
      as.character(observation_id),
      as.character(observation_id)
    )
  )

log_trace('Instantiating {round(as.double(nrow(observations)) * nrow(basis_vectors) * 8 / 1024 / 1024)} MB matrix')
H <- matrix(0, nrow(observations), nrow(basis_vectors))

for (sensitivities_base in args$sensitivities_base) {
  log_trace('Loading sensitivities from {sensitivities_base}')
  process_sensitivities(
    args$control,
    sensitivities_base,
    basis_vectors,
    levels(observations$observation_id),
    function(sensitivities, is_done) {
      log_trace('Filling matrix')
      H[cbind(
        to_integer(sensitivities$observation_id[!is_done]),
        to_integer(sensitivities$basis_vector[!is_done])
      )] <<- sensitivities$value[!is_done]
    }
  )
}

log_info('Writing output to {args$output}')
writeBin(
  as.vector(H),
  pipe(sprintf('lz4 -f - %s', args$output), 'wb')
)

log_info('Done')
