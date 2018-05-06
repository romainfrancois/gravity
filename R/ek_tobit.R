#' @title Eaton and Kortum (2001) Tobit model (EK Tobit)
#' 
#' @description \code{ek_tobit} estimates gravity models in their additive form
#' by conducting a censored regression.
#' It follows the Eaton and Kortum (2001) Tobit model where each country 
#' is assigned specific ceonsoring bounds.
#' 
#' @details \code{ek_tobit} represents the Eaton and Kortum (2001) Tobit model.
#' When taking the log of the gravity equation flows equal to zero
#' constitute a problem as their log is not defined.
#' Therefore, in \code{ek_tobit} all values of the dependent variable 
#' are redefined as intervals.
#' The positive observations have both interval bounds equal
#' to their original value. 
#' For zero flows the interval is left open. The right border
#' of the interval is set to the log of the minimum positive trade flow of 
#' the respective importing country.
#' The defined data object of class \code{\link[survival]{Surv}} is
#' then inserted in \code{\link[survival]{survreg}} for the 
#' parameter estimation.
#'  
#' To execute the function a square gravity dataset with all pairs of 
#' countries, ISO-codes for the country of origin and destination, a measure of 
#' distance between the bilateral partners as well as all 
#' information that should be considered as dependent an independent 
#' variables is needed. 
#' Missing bilateral flows as well as incomplete rows should be 
#' excluded from the dataset.  
#' Zero trade flows are allowed. 
#' 
#' \code{ek_tobit} is designed to be consistent with the Stata code provided at
#' the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' 
#' Up to now, the function is designed for cross-sectional data,
#' but can be extended to panel data using the 
#' \code{\link[survival]{survreg}} function.
#' 
#' For other Tobit functions, see \code{\link[gravity]{tobit}}
#' for a simple Tobit model where number \code{1} is added to all observations
#' and \code{\link[gravity]{et_tobit}} for the Eaton and Tamura (1994) 
#' threshold Tobit model where instead of simply adding number \code{1} 
#' to the data the threshold is estimated.
#' 
#' @param dependent_variable name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows.
#' The variable is logged and then taken as the dependent variable in 
#' the regression. As the log of zero is not defined, 
#' all flows equal to zero are replaced by
#' a left open interval with the logged minimum trade flow of the
#' respective importing country as right border.
#' 
#' @param regressors name (type: character) of the distance variable in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners. It is logged automatically when the function is executed. 
#' Interaction terms can be added.
#' 
#' @param code_destination variable name (type: character) of the label of the country 
#' (i.e ISO code) of destination in the dataset \code{data}. The variables 
#' are grouped by using \code{code_d} to obtain estimates.
#' 
#' @param vce_robust robust (type: logic) determines whether a robust 
#' variance-covariance matrix should be used. The default is set to \code{TRUE}. 
#' If set \code{TRUE} the estimation results are consistent with the 
#' Stata code provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' 
#' @param data name of the dataset to be used (type: character). 
#' To estimate gravity equations, a square gravity dataset including bilateral 
#' flows defined by the argument \code{y}, ISO-codes of type character 
#' (called \code{iso_o} for the country of origin and \code{iso_d} for the 
#' destination country), a distance measure defined by the argument \code{dist} 
#' and other potential influences given as a vector in \code{x} are required. 
#' All dummy variables should be of type numeric (0/1). Missing trade flows as 
#' well as incomplete rows should be excluded from the dataset. 
#' Zero trade flows are allowed.
#' 
#' @param ... additional arguments to be passed to \code{ek_tobit}.
#' 
#' @references 
#' For more information on gravity models, theoretical foundations and
#' estimation methods in general see 
#' 
#' Anderson, J. E. (1979) <DOI:10.12691/wjssh-2-2-5>
#' 
#' Anderson, J. E. (2010) <DOI:10.3386/w16576>
#' 
#' Anderson, J. E. and van Wincoop, E. (2003) <DOI:10.3386/w8079> 
#' 
#' Baier, S. L. and Bergstrand, J. H. (2009) <DOI:10.1016/j.jinteco.2008.10.004>
#' 
#' Baier, S. L. and Bergstrand, J. H. (2010) in Van Bergeijk, P. A., & Brakman, S. (Eds.) (2010) chapter 4 <DOI:10.1111/j.1467-9396.2011.01000.x>
#' 
#' Head, K., Mayer, T., & Ries, J. (2010) <DOI:10.1016/j.jinteco.2010.01.002>
#' 
#' Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>
#' 
#' Santos-Silva, J. M. C. and Tenreyro, S. (2006) <DOI:10.1162/rest.88.4.641> 
#' 
#' and the citations therein.
#' 
#' 
#' Especially for Tobit models see
#' 
#' Tobin, J. (1958) <DOI:10.2307/1907382>
#' 
#' Eaton, J., & Tamura, A. (1994) <DOI:10.3386/w4758>
#' 
#' Eaton, J., & Kortum, S. (2001) <DOI:10.3386/w8070>.
#' 
#' 
#' See Carson, R. T., & Sun, Yixiao (2007) <DOI:10.1111/j.1368-423X.2007.00218.x>
#' for the estimation of the threshold in a Tobit model. 
#' 
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' @examples 
#' \dontrun{
#' # Example for data with zero trade flows
#' data(gravity_zeros)
#' 
#' gravity_zeros <- gravity_zeros %>%
#'     mutate(
#'         lgdp_o = log(gdp_o),
#'         lgdp_d = log(gdp_d)
#'     )
#' 
#' ek_tobit(dependent_variable = "flow", regressors = c("distw", "rta","lgdp_o","lgdp_d"), 
#' code_destination = "iso_d",
#' vce_robust = TRUE, data = gravity_zeros)
#' }
#' 
#' \dontshow{
#' # examples for CRAN checks:
#' # executable in < 5 sec together with the examples above
#' # not shown to users
#' 
#' data(gravity_zeros)
#' gravity_zeros$lgdp_o <- log(gravity_zeros$gdp_o)
#' gravity_zeros$lgdp_d <- log(gravity_zeros$gdp_d)
#' 
#' # choose exemplarily 10 biggest countries for check data
#' countries_chosen_zeros <- names(sort(table(gravity_zeros$iso_o), decreasing = TRUE)[1:10])
#' grav_small_zeros <- gravity_zeros[gravity_zeros$iso_o %in% countries_chosen_zeros,]
#' ek_tobit(dependent_variable = "flow", regressors = c("distw", "rta", "lgdp_o", "lgdp_d"), 
#' code_destination = "iso_d",
#' vce_robust = TRUE, data = grav_small_zeros)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as a 
#' \code{\link[survival]{survreg}}-object.
#' 
#' @seealso \code{\link[survival]{Surv}}, \code{\link[survival]{survreg}}
#' 
#' @export 

ek_tobit <- function(dependent_variable, regressors, code_destination, vce_robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(vce_robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(code_destination) | code_destination %in% colnames(data) | length(code_destination) == 1)
  
  # Split input vectors -----------------------------------------------------
  distance <- regressors[1]
  additional_regressors <- regressors[-1]
  
  # Discarding unusable observations ----------------------------------------
  d <- data %>% 
    filter_at(vars(!!sym(distance)), any_vars(!!sym(distance) > 0)) %>% 
    filter_at(vars(!!sym(distance)), any_vars(is.finite(!!sym(distance))))
  
  # Transforming data, logging distances ---------------------------------------
  d <- d %>% 
    mutate(
      dist_log = log(!!sym(distance))
    )
  
  # Transforming data, logging flows -------------------------------------------
  d <- d %>% 
    mutate(
      y_log_ek = ifelse(!!sym(dependent_variable) > 0, log(!!sym(dependent_variable)), NA)
    )
  
  # Minimum flows -----------------------------------------------------------
  d <- d %>% 
    group_by(!!sym(code_destination)) %>% 
    mutate(
      exportmin = min(!!sym("y_log_ek"), na.rm = TRUE)
    )
  
  # Transforming censored variables -----------------------------------------
  d <- d %>% 
    mutate(
      flows_ek1 = ifelse(!!sym(dependent_variable) > 0, !!sym("y_log_ek"), -Inf)
    ) %>% 
    
    mutate(
      flows_ek2 = ifelse(!!sym(dependent_variable) > 0, !!sym("flows_ek1"), !!sym("exportmin"))
    ) %>% 
    
    ungroup()
  
  # Response variable -------------------------------------------------------
  f1 <- d %>% select(!!sym("flows_ek1")) %>% as_vector()
  f2 <- d %>% select(!!sym("flows_ek2")) %>% as_vector()
  
  y_cens_log_ek <- survival::Surv(f1, f2, type = "interval2") %>% as_vector()
  
  # Model -------------------------------------------------------------------
  vars           <- paste(c("dist_log", additional_regressors), collapse = " + ")
  form           <- stats::as.formula(paste("y_cens_log_ek", "~", vars))
  model_ek_tobit <- survival::survreg(form, data = d, dist = "gaussian", robust = vce_robust)
  
  # Return ------------------------------------------------------------------
  return_object      <- summary(model_ek_tobit)
  return_object$call <- form
  return(return_object)
}
