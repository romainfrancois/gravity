#' @title Double Demeaning (DDM)
#' 
#' @description \code{ddm} estimates gravity models via double demeaning the 
#' left hand side and right hand side of the gravity equation.
#' 
#' @details \code{ddm} is an estimation method for gravity models presented
#' in Head and Mayer (2014) (see the references for more information).  
#' To execute the function a square gravity dataset with all pairs of 
#' countries, ISO-codes for the country of origin and destination, a measure of 
#' distance between the bilateral partners as well as all 
#' information that should be considered as dependent an independent 
#' variables is needed. 
#' Make sure the ISO-codes are of type "character".
#' Missing bilateral flows as well as incomplete rows should be 
#' excluded from the dataset. 
#' Furthermore, flows equal to zero should be excluded as the gravity equation 
#' is estimated in its additive form.  
#' Country specific effects are subdued due double demeaning. 
#' Hence, unilateral income proxies such as GDP cannot be 
#' considered as exogenous variables.
#' 
#' \code{ddm} is designed to be consistent with the Stata code provided at
#' the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' As, to our knowledge at the moment, there is no explicit literature covering 
#' the estimation of a gravity equation by \code{ddm} using panel data, 
#' we do not recommend to apply this method in this case.
#' 
#' @param dependent_variable name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. This variable is logged and then used as the 
#' dependent variable in the estimation.
#' 
#' @param regressors name (type: character) of the distance variable in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners. It is logged automatically when the function is executed. 
#' 
#' @param codes variable name (type: character) of the label of the country 
#' (i.e ISO code) of origin in the dataset \code{data}. The variables 
#' are grouped by using \code{code_o} and \code{code_d} to obtain estimates. 
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
#' Furthermore, flows equal to zero should be excluded as the gravity equation 
#' is estimated in its additive form.
#' As, to our knowledge at the moment, there is no explicit literature covering 
#' the estimation of a gravity equation by \code{ddm} 
#' using panel data, cross-sectional data should be used. 
#' 
#' @param ... additional arguments to be passed to \code{ddm}.
#' 
#' @references 
#' For more information on Double Demeaning as well as information on gravity 
#' models, theoretical foundations and estimation methods in general see
#' 
#' Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>
#' 
#' as well as
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
#' Santos-Silva, J. M. C. and Tenreyro, S. (2006) <DOI:10.1162/rest.88.4.641> 
#' 
#' and the citations therein.
#' 
#' 
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' @examples 
#' \dontrun{
#' data(gravity_no_zeros)
#' 
#' ddm(dependent_variable = "flow", regressors = c("distw", "rta"), 
#' codes = c("iso_o", "iso_d"), 
#' vce_robust = TRUE, data = gravity_no_zeros)
#' 
#' ddm(dependent_variable = "flow", regressors = c("distw", "rta", "comcur", "contig"), 
#' codes = c("iso_o", "iso_d"), 
#' vce_robust=TRUE, data=gravity_no_zeros)
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
#' ddm(dependent_variable = "flow", regressors = c("distw", "rta"), 
#' codes = c("iso_o", "iso_d"), 
#' vce_robust = TRUE, data = grav_small)
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

ddm <- function(dependent_variable, regressors, codes, vce_robust=TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(vce_robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(codes) | all(codes %in% colnames(data)) | length(codes) == 2)
  
  # Split input vectors -----------------------------------------------------
  code_o <- codes[1]
  code_d <- codes[2]
  
  distance <- regressors[1]
  additional_regressors <- regressors[-1]
  
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
    
    group_by(!!sym(code_o), add = FALSE) %>% 
    mutate(
      ym1 = mean(!!sym("y_log_ddm"), na.rm = TRUE),
      dm1 = mean(!!sym("dist_log_ddm"), na.rm = TRUE)
    ) %>% 
    
    group_by(!!sym(code_d), add = FALSE) %>% 
    mutate(
      ym2 = mean(!!sym("y_log_ddm"), na.rm = TRUE),
      dm2 = mean(!!sym("dist_log_ddm"), na.rm = TRUE)
    ) %>% 
    
    group_by(!!sym(code_o), add = FALSE) %>% 
    mutate(
      y_log_ddm = !!sym("y_log_ddm") - !!sym("ym1"),
      dist_log_ddm = !!sym("dist_log_ddm") - !!sym("dm1")
    ) %>% 
    
    group_by(!!sym(code_d), add = FALSE) %>% 
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
    select(!!sym(code_o), !!sym(code_d), additional_regressors) %>% 
    gather(!!sym("key"), !!sym("value"), -!!sym(code_o), -!!sym(code_d)) %>% 
    
    mutate(key = paste0(!!sym("key"), "_ddm")) %>% 
    
    group_by(!!sym(code_o), !!sym("key"), add = FALSE) %>% 
    mutate(ddm = !!sym("value") - mean(!!sym("value"), na.rm = TRUE)) %>% 
    
    group_by(!!sym(code_d), !!sym("key"), add = FALSE) %>% 
    mutate(ddm = !!sym("ddm") - mean(!!sym("value"), na.rm = TRUE)) %>% 
    
    ungroup() %>% 
    mutate(value = !!sym("ddm") + mean(!!sym("value"), na.rm = TRUE)) %>% 
    
    select(!!!syms(c(code_o, code_d, "key", "value"))) %>% 
    spread(!!sym("key"), !!sym("value"))
  
  # Model ----------------------------------------------------------------------
  dmodel <- left_join(d, d2, by = c(code_o, code_d)) %>% 
    select(!!sym("y_log_ddm"), ends_with("_ddm"))
  
  model_ddm <- stats::lm(y_log_ddm ~ . + 0, data = dmodel)
  
  # Return ---------------------------------------------------------------------
  if (vce_robust == TRUE) {
    return_object      <- robust_summary(model_ddm, robust = TRUE)
    return_object$call <- as.formula(model_ddm)
    return(return_object)
  }
  
  if (vce_robust == FALSE) {
    return_object      <- robust_summary(model_ddm, robust = FALSE)
    return_object$call <- as.formula(model_ddm)
    return(return_object)
  }
}
