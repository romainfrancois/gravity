#' @title Double Demeaning (DDM)
#'
#' @description \code{ddm} estimates gravity models via double demeaning the
#' left hand side and right hand side of the gravity equation.
#'
#' @details \code{ddm} is an estimation method for gravity models presented
#' in \insertCite{Head2014;textual}{gravity}.
#'
#' Country specific effects are subdued due double demeaning. Hence, unilateral income
#' proxies such as GDP cannot be considered as exogenous variables.
#'
#' Unilateral effect drop out due to double demeaning and therefore cannot be estimated.
#'
#' \code{ddm} is designed to be consistent with the Stata code provided at
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#'
#' As, to our knowledge at the moment, there is no explicit literature covering
#' the estimation of a gravity equation by \code{ddm} using panel data,
#' we do not recommend to apply this method in this case.
#'
#' @param dependent_variable (Type: character) name of the dependent variable in the dataset
#' \code{data} (e.g. trade flows).
#'
#' This dependent variable is divided by the
#' product of unilateral incomes (e.g.
#' GDPs or GNPs of the countries of interest, named \code{inc_o} and \code{inc_d} in the example datasets)
#' and logged afterwards.
#'
#' The transformed variable is then used as the dependent variable in the
#' estimation.
#'
#' @param distance (Type: character) name of the distance variable in the dataset \code{data} containing a measure of
#' distance between all pairs of bilateral partners and bilateral variables that should
#' be taken as the independent variables in the estimation.
#'
#' The distance is logged automatically when the function is executed.
#'
#' @param additional_regressors (Type: character) names of the additional regressors to include in the model (e.g. a dummy 
#' variable to indicate contiguity).
#'
#' Unilateral metric variables such as GDPs should be inserted via the arguments \code{income_origin} and \code{income_origin}.
#' 
#' As country specific effects are subdued due to demeaning, no further unilateral variables apart from incomes can be added.
#'
#' Write this argument as \code{c(contiguity, common currency, ...)}.
#'
#' @param code_origin (Type: character) variable name of the code of the country
#' of origin (e.g. ISO-3 codes from the variables \code{iso_o} in the
#' example datasets). The variables are grouped by using \code{iso_o} and \code{iso_d} to obtain estimates.
#'
#' @param code_destination (Type: character) variable name of the code of the country
#' of destination (e.g. ISO-3 codes from the variables \code{iso_d}) in the
#' example datasets). The variables are grouped by using \code{iso_o} and \code{iso_d} to obtain estimates.
#'
#' @param robust (Type: logical) determines whether robust
#' fitting should be used. By default is set to \code{TRUE}.
#' 
#' @param data (Type: character) name of the dataset to be used.
#'
#' To estimate gravity equations you need a square dataset including bilateral
#' flows defined by the argument \code{dependent_variable}, ISO codes or similar of type character
#' (e.g. \code{iso_o} for the country of origin and \code{iso_d} for the
#' destination country), a distance measure defined by the argument \code{distance}
#' and other potential influences (e.g. contiguity and common currency) given as a vector in
#' \code{regressors} are required.
#'
#' All dummy variables should be of type numeric (0/1).
#'
#' Make sure the ISO codes are of type "character".
#'
#' If an independent variable is defined as a ratio, it should be logged.
#'
#' The user should perform some data cleaning beforehand to remove observations that contain entries that
#' can distort estimates.
#'
#' The function will remove zero flows and distances.
#'
#' @param ... additional arguments to be passed to \code{ddm}.
#'
#' @references
#' For more information on Double Demeaning as well as information on gravity
#' models, theoretical foundations and estimation methods in general see
#'
#' \insertRef{Head2014}{gravity}
#'
#' as well as
#'
#' \insertRef{Anderson1979}{gravity}
#'
#' \insertRef{Anderson2001}{gravity}
#'
#' \insertRef{Anderson2010}{gravity}
#'
#' \insertRef{Baier2009}{gravity}
#'
#' \insertRef{Baier2010}{gravity}
#'
#' \insertRef{Head2010}{gravity}
#'
#' \insertRef{Santos2006}{gravity}
#'
#' and the citations therein.
#'
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} 
#' for gravity datasets and Stata code for estimating gravity models.
#'
#' @examples
#' \dontrun{
#' data(gravity_no_zeros)
#'
#' ddm(dependent_variable = "flow", 
#' distance = "distw", additional_regressors = "rta",
#' code_origin = "iso_o", code_destination = "iso_d",
#' robust = TRUE, data = gravity_no_zeros)
#'
#' ddm(dependent_variable = "flow", 
#' distance = "distw", additional_regressors = c("distw", "rta", "comcur", "contig"),
#' code_origin = "iso_o", code_destination = "iso_d",
#' robust = TRUE, data = gravity_no_zeros)
#' }
#'
#' \dontshow{
#' # examples for CRAN checks:
#' # executable in < 5 sec together with the examples above
#' # not shown to users
#'
#' data(gravity_no_zeros)
#' # choose exemplarily 10 biggest countries for check data
#' countries_chosen <- names(sort(table(gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
#' grav_small <- gravity_no_zeros[gravity_no_zeros$iso_o %in% countries_chosen,]
#' ddm(dependent_variable = "flow", distance = "distw", additional_regressors = "rta",
#' code_origin = "iso_o", code_destination = "iso_d",
#' robust = TRUE, data = grav_small)
#' }
#'
#' @return
#' The function returns the summary of the estimated gravity model as an
#' \code{\link[stats]{lm}}-object.
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}},
#' \code{\link[sandwich]{vcovHC}}
#'
#' @export

ddm <- function(dependent_variable = "flow", distance = "distw", additional_regressors = NULL, 
                code_origin = "iso_o", code_destination = "iso_d", 
                robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(robust))
  
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  
  stopifnot(is.character(distance), distance %in% colnames(data), length(distance) == 1)
  
  if (!is.null(additional_regressors)) {
    stopifnot(is.character(additional_regressors), all(additional_regressors %in% colnames(data)))
  }
  
  stopifnot(is.character(code_origin) | code_origin %in% colnames(data) | length(code_origin) == 1)
  stopifnot(is.character(code_destination) | code_destination %in% colnames(data) | length(code_destination) == 1)
  
  # Discarding unusable observations ----------------------------------------
  d <- data %>%
    filter_at(vars(!!sym(distance)), any_vars(!!sym(distance) > 0)) %>%
    filter_at(vars(!!sym(distance)), any_vars(is.finite(!!sym(distance)))) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(!!sym(dependent_variable) > 0)) %>%
    filter_at(vars(!!sym(dependent_variable)), any_vars(is.finite(!!sym(dependent_variable))))

  # Transforming data, logging distances ---------------------------------------
  d <- d %>%
    mutate(
      dist_log = log(!!sym(distance))
    )

  # Transforming data, logging flows -------------------------------------------
  d <- d %>%
    mutate(
      y_log = log(!!sym(dependent_variable))
    )

  # Substracting the means -----------------------------------------------------
  d <- d %>%
    mutate(
      y_log_ddm = !!sym("y_log"),
      dist_log_ddm = !!sym("dist_log")
    ) %>%
    group_by(!!sym(code_origin), add = FALSE) %>%
    mutate(
      ym1 = mean(!!sym("y_log_ddm"), na.rm = TRUE),
      dm1 = mean(!!sym("dist_log_ddm"), na.rm = TRUE)
    ) %>%
    group_by(!!sym(code_destination), add = FALSE) %>%
    mutate(
      ym2 = mean(!!sym("y_log_ddm"), na.rm = TRUE),
      dm2 = mean(!!sym("dist_log_ddm"), na.rm = TRUE)
    ) %>%
    group_by(!!sym(code_origin), add = FALSE) %>%
    mutate(
      y_log_ddm = !!sym("y_log_ddm") - !!sym("ym1"),
      dist_log_ddm = !!sym("dist_log_ddm") - !!sym("dm1")
    ) %>%
    group_by(!!sym(code_destination), add = FALSE) %>%
    mutate(
      y_log_ddm = !!sym("y_log_ddm") - !!sym("ym2"),
      dist_log_ddm = !!sym("dist_log_ddm") - !!sym("dm2")
    ) %>%
    ungroup() %>%
    mutate(
      y_log_ddm = !!sym("y_log_ddm") + mean(!!sym("y_log"), na.rm = TRUE),
      dist_log_ddm = !!sym("dist_log_ddm") + mean(!!sym("dist_log"), na.rm = TRUE)
    )

  # Substracting the means for the other independent variables -----------------
  d2 <- d %>%
    select(!!sym(code_origin), !!sym(code_destination), !!!syms(additional_regressors)) %>%
    gather(!!sym("key"), !!sym("value"), -!!sym(code_origin), -!!sym(code_destination)) %>%
    mutate(key = paste0(!!sym("key"), "_ddm")) %>%
    group_by(!!sym(code_origin), !!sym("key"), add = FALSE) %>%
    mutate(ddm = !!sym("value") - mean(!!sym("value"), na.rm = TRUE)) %>%
    group_by(!!sym(code_destination), !!sym("key"), add = FALSE) %>%
    mutate(ddm = !!sym("ddm") - mean(!!sym("value"), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(value = !!sym("ddm") + mean(!!sym("value"), na.rm = TRUE)) %>%
    select(!!!syms(c(code_origin, code_destination, "key", "value"))) %>%
    spread(!!sym("key"), !!sym("value"))

  # Model ----------------------------------------------------------------------
  if (!is.null(additional_regressors)) {
    d <- left_join(d, d2, by = c(code_origin, code_destination)) %>%
      select(!!sym("y_log_ddm"), ends_with("_ddm"))
  } else {
    d <- select(d, !!sym("y_log_ddm,"), ends_with("_ddm"))
  }
  
  if (robust == TRUE) {
    model_bvw <- MASS::rlm(y_log_ddm ~ . + 0, data = d)
  } else {
    model_bvw <- stats::lm(y_log_ddm ~ . + 0, data = d)
  }
  
  return(model_bvw)
}
