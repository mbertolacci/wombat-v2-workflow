library(argparse)
library(dplyr, warn.conflicts = FALSE)

source(Sys.getenv('UTILS_PARTIAL'))

to_basis_vector_str <- function(x) {
  with(x, {
    sprintf(
      '%s.%s.%s%s%s',
      as.character(inventory),
      as.character(component),
      as.character(region),
      ifelse(is.na(month), '', '.'),
      ifelse(is.na(month), '', as.character(month))
    )
  })
}

parser <- ArgumentParser()
parser$add_argument('--perturbations')
parser$add_argument('--output')
args <- parser$parse_args()

perturbations <- fst::read_fst(args$perturbations)
basis_vector_df <- bind_rows(
  expand.grid(
    region = levels(perturbations$region),
    inventory = levels(perturbations$inventory),
    component = levels(perturbations$component)
  ) %>%
    filter(
      component != 'residual',
      !(inventory == 'ocean' & grepl('3', component))
    ),
  expand.grid(
    region = levels(perturbations$region),
    inventory = levels(perturbations$inventory),
    month = levels(perturbations$month)
  ) %>%
    mutate(component = factor('residual'))
) %>%
  mutate(
    basis_vector_str = to_basis_vector_str(.),
    basis_vector = factor(basis_vector_str, basis_vector_str)
  )

fst::write_fst(basis_vector_df, args$output)
