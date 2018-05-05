#' @title Bonus vetus OLS (BVU)
#' 
#' @description \code{bvu} estimates gravity models via Bonus vetus OLS with simple averages.
#' 
#' @details \code{Bonus vetus OLS} is an estimation method for gravity models 
#' developed by Baier and Bergstrand (2009, 2010) using simple averages to center a
#' Taylor-series (see the references for more information). 
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
#' The \code{bvu} function considers Multilateral Resistance terms and allows to 
#' conduct comparative statics. Country specific effects are subdued due 
#' to demeaning. Hence, unilateral variables apart from \code{inc_o}
#' and \code{inc_d} cannot be included in the estimation.
#' \code{bvu} is designed to be consistent with the Stata code provided at
#' the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' As, to our knowledge at the moment, there is no explicit literature covering 
#' the estimation of a gravity equation by \code{bvu} using panel data, 
#' we do not recommend to apply this method in this case.
#' 
#' @param dependent_variable name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. This dependent variable is divided by the 
#' product of unilateral incomes (named \code{inc_o} and \code{inc_d}, e.g. 
#' GDPs or GNPs of the countries of interest) and logged afterwards.
#' The transformed variable is then used as the dependent variable in the 
#' estimation.
#' 
#' @param regressors name (type: character) of the distance variable () in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners and bilateral variables that should be taken as the independent variables 
#' in the estimation. The distance is logged automatically when the function is executed. 
#' If an independent variable is a dummy variable, it should be of type numeric (0/1) in the dataset. 
#' If an independent variable is defined as a ratio, it should be logged. 
#' Unilateral metric variables such as GDPs should be inserted via the argument \code{incomes} 
#' As country specific effects are subdued due to demeaning, no further unilateral variables 
#' equivalent to \code{c(gdp origin, gdp destination)} can be added.
#' Write this argument as \code{c(distance, other variables)}. 
#' 
#' @param incomes variable name (type: character) of the income of the country of 
#' origin (i.e \code{gdp_o}) and destination (i.e \code{gdp_d}) in the dataset \code{data}. 
#' The dependent variable \code{dependent_variable} is divided by the product of the incomes.
#' Write this argument as \code{c(origin, destination)}. 
#' 
#' @param codes variable name (type: character) of the code of the country 
#' (i.e ISO-3 code) of origin (i.e \code{iso_o}) and destination (i.e \code{iso_d}) in the dataset \code{data}. 
#' The variables are grouped by using \code{code_o} and \code{code_d} to obtain estimates.
#' Write this argument as \code{c(origin, destination)}.
#' 
#' @param vce_robust robust (type: logical) determines whether a robust 
#' variance-covariance matrix should be used. The default is set to \code{TRUE}. 
#' If set \code{TRUE} the estimation results are consistent with the 
#' Stata code provided at the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' 
#' @param data name of the dataset to be used (type: character). 
#' As an example, to estimate gravity equations you need a square dataset including bilateral 
#' flows defined by the argument \code{dependent_variable}, ISO-codes of type character 
#' (i.e. \code{iso_o} for the country of origin and \code{iso_d} for the 
#' destination country), a distance measure defined by the argument \code{distance} 
#' and other potential influences (i.e. contiguity or common currency) given as a vector in 
#' \code{regressors} are required. 
#' All dummy variables should be of type numeric (0/1). Some of our functions remove observations with 
#' null/missing trade flows or distances as those variables are converted to log scale before obtaining estimates. 
#' The user should perform some data cleaning beforehand to remove observations that contain entries that 
#' can distort estimates.
#' As, to our knowledge at the moment, there is no explicit literature covering the estimation of a gravity 
#' equation by \code{bvu} using panel data, cross-sectional data should be used. 
#' 
#' @param ... additional arguments to be passed to \code{bvu}.
#' 
#' @references 
#' For estimating gravity equations via Bonus Vetus OLS see
#' 
#' Baier, S. L. and Bergstrand, J. H. (2009) <DOI:10.1016/j.jinteco.2008.10.004>
#' 
#' Baier, S. L. and Bergstrand, J. H. (2010) in Van Bergeijk, P. A., & Brakman, S. (Eds.) (2010) chapter 4 <DOI:10.1111/j.1467-9396.2011.01000.x>
#' 
#' For more information on gravity models, theoretical foundations and
#' estimation methods in general see 
#' 
#' Anderson, J. E. (1979) <DOI:10.12691/wjssh-2-2-5>
#' 
#' Anderson, J. E. (2010) <DOI:10.3386/w16576>
#' 
#' Anderson, J. E. and van Wincoop, E. (2003) <DOI:10.3386/w8079> 
#' 
#' Head, K., Mayer, T., & Ries, J. (2010) <DOI:10.1016/j.jinteco.2010.01.002>
#' 
#' Head, K. and Mayer, T. (2014) <DOI:10.1016/B978-0-444-54314-1.00003-3>
#' 
#' Santos-Silva, J. M. C. and Tenreyro, S. (2006) <DOI:10.1162/rest.88.4.641> 
#' 
#' and the citations therein.
#' 
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' @examples 
#' \dontrun{
#' data(gravity_no_zeros)
#' 
#' bvu(dependent_variable = "flow", regressors = c("distw", "rta"), 
#' incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
#' vce_robust = TRUE, data = gravity_no_zeros)
#' 
#' bvu(dependent_variable = "flow", regressors = c("distw", "rta", "contig", "comcur"), 
#' incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
#' vce_robust = TRUE, data = gravity_no_zeros)
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
#' bvu(dependent_variable = "flow", regressors = c("distw", "rta"), 
#' incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
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

bvu <- function(dependent_variable, regressors, incomes, codes, vce_robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(vce_robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)))
  stopifnot(is.character(incomes) | all(incomes %in% colnames(data)) | length(incomes) == 2)
  stopifnot(is.character(codes) | all(codes %in% colnames(data)) | length(codes) == 2)

  # Split input vectors -----------------------------------------------------
  inc_o <- incomes[[1]]
  inc_d <- incomes[[2]]
  
  code_o <- codes[[1]]
  code_d <- codes[[2]]
  
  distance <- regressors[[1]]
  additional_regressors <- regressors[[-1]]
  
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
      y_log_bvu = log(
        !!sym(dependent_variable) / (!!sym(inc_o) * !!sym(inc_d))
      )
    )
  
  # Multilateral Resistance (MR) for distance ----------------------------------
  d <- d %>% 
    group_by(!!sym(code_o)) %>% 
    mutate(mean_dist_log_1 = mean(!!sym("dist_log"), na.rm = TRUE)) %>% 
    group_by(!!sym(code_d), add = FALSE) %>% 
    mutate(mean_dist_log_2 = mean(!!sym("dist_log"), na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(
      mean_dist_log_3 = mean(!!sym("dist_log"), na.rm = TRUE),
      dist_log_mr = !!sym("dist_log") - 
        (!!sym("mean_dist_log_1") + !!sym("mean_dist_log_2") - !!sym("mean_dist_log_3"))
    )
  
  # Multilateral Resistance (MR) for the other independent variables -----------
  d2 <- d %>% 
    select(!!sym(code_o), !!sym(code_d), additional_regressors) %>% 
    gather(!!sym("key"), !!sym("value"), -!!sym(code_o), -!!sym(code_d)) %>% 
    
    group_by(!!sym(code_o), !!sym("key")) %>% 
    mutate(mean_dist_log_1 = mean(!!sym("value"), na.rm = TRUE)) %>% 
    
    group_by(!!sym(code_d), !!sym("key")) %>% 
    mutate(mean_dist_log_2 = mean(!!sym("value"), na.rm = TRUE)) %>% 
    
    group_by(!!sym("key")) %>% 
    mutate(
      mean_dist_log_3 = mean(!!sym("value"), na.rm = TRUE),
      dist_log_mr = !!sym("value") - (!!sym("mean_dist_log_1") + !!sym("mean_dist_log_2") - !!sym("mean_dist_log_3"))
    ) %>% 
    
    ungroup() %>% 
    mutate(key = paste0(!!sym("key"), "_mr")) %>% 
    
    select(!!!syms(c(code_o, code_d, "key", "dist_log_mr"))) %>% 
    spread(!!sym("key"), !!sym("dist_log_mr"))
  
  # Model ----------------------------------------------------------------------
  dmodel <- left_join(d, d2, by = c(code_o, code_d)) %>% 
    select(!!sym("y_log_bvu"), ends_with("_mr"))
  
  model_bvu <- stats::lm(y_log_bvu ~ ., data = dmodel)
  
  # Return ---------------------------------------------------------------------
  if (vce_robust == TRUE) {
    return_object      <- robust_summary(model_bvu, robust = TRUE)
    return_object$call <- as.formula(model_bvu)
    return(return_object)
  }
  
  if (vce_robust == FALSE) {
    return_object      <- robust_summary(model_bvu, robust = FALSE)
    return_object$call <- as.formula(model_bvu)
    return(return_object)
  }
}
