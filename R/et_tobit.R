#' @title Eaton and Tamura (1994) threshold Tobit model (ET Tobit)
#' 
#' @description \code{et_tobit} estimates gravity models in their additive form
#' by conducting a left-censored regression.
#' It follows the Eaton and Tamura (1994) Tobit model,
#' also called threshold Tobit model, where,
#' instead of adding number \code{1} to the dependent variable as done 
#' in \code{\link[gravity]{tobit}}, the constant added to the
#' data is estimated and interpreted as a threshold.
#' For estimating this threshold, we follow Carson and Sun (2007).
#' 
#' @details \code{et_tobit} represents the Eaton and Tamura (1994) Tobit model
#' which is often used when several gravity models are compared.
#' When taking the log of the gravity equation flows equal to zero
#' constitute a problem as their log is not defined.
#' Therefore, a constant is added to the flows. 
#' This constant, opposed to \code{\link[gravity]{tobit}}, is estimated. 
#' Compared to the usual ET-Tobit approaches, in this package, the estimation
#' of the threshold is done before the other parameters are estimated.
#' We follow Carson and Sun (2007), who show that taking the minimum
#' positive flow value as an estimate of the threshold is super-consistent and that
#' using this threshold estimate ensures that the parameter MLE 
#' are asymptotically normal with the asymptotic variance
#' identical to the variance achieved when the threshold is known. 
#' Hence, first the threshold is estimated as the minimum positive flow. 
#' This threshold is added to the flow variable. It is logged
#' afterwards and taken as the dependent variable. 
#' The Tobit estimation is then conducted using the 
#' \code{\link[censReg]{censReg}} function and setting the lower bound 
#' equal to the log of the minimum positive flow value which was added to all
#' observations.
#' A Tobit regression represents a combination of a binary and a
#' linear regression. 
#' This procedure has to be taken into consideration when
#' interpreting the estimated coefficients.
#' The marginal effects of an explanatory variable on the expected value of 
#' the dependent variable equals the product of both the probability of the 
#' latent variable exceeding the threshold and the marginal effect of the 
#' explanatory variable of the expected value of the latent variable. 
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
#' Up to now, the function is designed for cross-sectional data,
#' but can be easily extended to panel data using the 
#' \code{\link[censReg]{censReg}} function.
#' A robust estimations is not implemented to the present
#' as the \code{\link[censReg]{censReg}} function is not
#' compatible with the \code{\link[sandwich]{vcovHC}} function.
#' 
#' For a more elaborate Tobit function, see \code{\link[gravity]{ek_tobit}} 
#' for the
#' Eaton and Kortum (2001) Tobit model where each zero trade volume
#' is assigned a country specific interval with the upper 
#' bound equal to the minimum positive trade level of the respective
#' importing country.
#' 
#' @param dependent_variable name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows.
#' Following Carson and Sun (2007), the smallest positive flow value is 
#' used as an estimate of the threshold. 
#' It is added to \code{y}, the transformed variable is 
#' logged and taken as the dependent variable in the Tobit estimation with 
#' lower bound equal to the log of the smallest possible flow value. 
#' 
#' @param regressors name (type: character) of the distance variable in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners. It is logged automatically when the function is executed. 
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
#' @param ... additional arguments to be passed to \code{et_tobit}.
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
#' et_tobit(y = "flow", regressors = c("distw", "rta","lgdp_o","lgdp_d"),
#' data = gravity_zeros)
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
#' et_tobit(dependent_variable = "flow", regressors = c("distw", "rta","lgdp_o","lgdp_d"),
#' data = grav_small_zeros)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as a 
#' \code{\link[censReg]{censReg}}-object.
#' 
#' @seealso \code{\link[censReg]{censReg}}
#' 
#' @export 

et_tobit <- function(dependent_variable, regressors, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  
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
  flow_min_log <- filter_at(d, vars(!!sym(dependent_variable)), any_vars(!!sym(dependent_variable) > 0))
  
  d <- d %>% 
    mutate(
      y_log_et = ifelse(!!sym(dependent_variable) > 0, log(!!sym(dependent_variable)), NA)
    )
  
  # Transforming data, logging flows, distances --------------------------------
  d <- d %>% 
    mutate(
      y2 = ifelse(!!sym(dependent_variable) > 0, !!sym(dependent_variable), NA),
      y2_log = log(!!sym("y2"))
    )
  
  y2min <- min(d %>% select(!!sym("y2")), na.rm = TRUE)
  y2min_log <- log(y2min)

  d <- d %>% 
    rowwise() %>% 
    mutate(y_cens_log_et = log(sum(!!sym(dependent_variable), y2min, na.rm = TRUE))) %>% 
    ungroup()
  
  # Model -------------------------------------------------------------------
  vars           <- paste(c("dist_log", additional_regressors), collapse = " + ")
  form           <- stats::as.formula(paste("y_cens_log_et", "~", vars))
  model_et_tobit <- censReg::censReg(formula = form, 
                                     left = y2min_log, right = Inf, 
                                     data = d, start = NULL)
  
  # Return ------------------------------------------------------------------
  return_object      <- summary(model_et_tobit)
  return_object$call <- form
  return(return_object)
}
