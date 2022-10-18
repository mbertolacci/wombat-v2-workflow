library(argparse)
library(dplyr, warn.conflicts = FALSE)
library(Rcpp)

source(Sys.getenv('UTILS_PARTIAL'))
sourceCpp(Sys.getenv('UTILS_CPP_PARTIAL'))

parser <- ArgumentParser()
parser$add_argument('--basis-vectors')
parser$add_argument('--control-emissions')
parser$add_argument('--perturbations')
parser$add_argument('--output')
args <- parser$parse_args()

basis_vectors <- fst::read_fst(args$basis_vectors)
control_emissions <- fst::read_fst(args$control_emissions)
perturbations <- fst::read_fst(args$perturbations)
cell_area <- control_emissions %>%
  distinct(longitude, latitude, cell_height, area) %>%
  mutate(
    latitude_bottom = latitude - cell_height / 2
  ) %>%
  select(-cell_height)

perturbations_base <- perturbations %>%
  add_basis_vector(basis_vectors) %>%
  mutate(
    minor_component = factor(
      case_when(
        component == 'residual' ~ 'residual',
        component %in% c('intercept', 'trend') ~ 'linear',
        TRUE ~ 'periodic'
      ),
      c('linear', 'periodic', 'residual')
    ),
    inventory_minor_component_time = interaction(
      inventory,
      minor_component,
      time,
      drop = TRUE
    )
  ) %>%
  left_join(cell_area, by = c('longitude', 'latitude'))

fst::write_fst(perturbations_base, args$output)
