#' gravity exported operators and S3 methods
#' The following functions are imported and then re-exported
#' from the gravity package to avoid listing Depends of gravity.
#' @importFrom dplyr select mutate group_by ungroup row_number left_join
#'   ends_with vars filter_at any_vars
#' @importFrom tidyr gather spread
#' @importFrom rlang sym syms quo
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @name gravity-exports
NULL