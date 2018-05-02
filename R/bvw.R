#' @title Bonus vetus OLS (BVW)
#' 
#' @description \code{bvw} estimates gravity models via Bonus 
#' vetus OLS with GDP-weights.
#' 
#' @details \code{Bonus vetus OLS} is an estimation method for gravity models 
#' developed by Baier and Bergstrand (2009, 2010) using GDP-weights to center a
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
#' The \code{bvw} function considers Multilateral Resistance terms and allows to 
#' conduct comparative statics. Country specific effects are subdued due 
#' to demeaning. Hence, unilateral variables apart from \code{inc_o}
#' and \code{inc_d} cannot be included in the estimation.
#' \code{bvw} is designed to be consistent with the Stata code provided at
#' the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' As, to our knowledge at the moment, there is no explicit literature covering 
#' the estimation of a gravity equation by \code{bvw} using panel data, 
#' we do not recommend to apply this method in this case.
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. This dependent variable is divided by the 
#' product of unilateral incomes (named \code{inc_o} and \code{inc_d}, e.g. 
#' GDPs or GNPs of the countries of interest) and logged afterwards.
#' The transformed variable is then used as the dependent variable in the 
#' estimation.
#' 
#' @param dist name (type: character) of the distance variable in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners. It is logged automatically when the function is executed. 
#' 
#' @param x vector of names (type: character) of those bilateral variables in 
#' the dataset \code{data} that should be taken as the independent variables 
#' in the estimation. If an independent variable is a dummy variable,
#' it should be of type numeric (0/1) in the dataset. If an independent variable 
#' is defined as a ratio, it should be logged. Unilateral metric variables 
#' such as GDPs should be inserted via the arguments \code{inc_o} 
#' for the country of origin and \code{inc_d} for the country of destination.
#' As country specific effects are subdued due to demeaning, no further
#' unilateral variables apart from \code{inc_o} and \code{inc_d} can be 
#' added.
#' 
#' @param inc_o variable name (type: character) of the income of the country of 
#' origin in the dataset \code{data}. The dependent variable \code{y} is
#' divided by the product of the incomes \code{inc_d} and \code{inc_o}. 
#' 
#' @param inc_d variable name (type: character) of the income of the country of 
#' destination in the dataset \code{data}. The dependent variable \code{y} is
#' divided by the product of the incomes \code{inc_d} and \code{inc_o}. 
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
#' the estimation of a gravity equation by \code{bvw} 
#' using panel data, cross-sectional data should be used. 
#' 
#' @param ... additional arguments to be passed to \code{bvw}.
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
#' bvw(y="flow", dist="distw", x=c("rta"), 
#' inc_o="gdp_o", inc_d="gdp_d", vce_robust=TRUE, data=gravity_no_zeros)
#' 
#' bvw(y="flow", dist="distw", x=c("rta", "comcur", "contig"), 
#' inc_o="gdp_o", inc_d="gdp_d", vce_robust=TRUE, data=gravity_no_zeros)
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
#' bvw(y="flow", dist="distw", x=c("rta"), inc_o="gdp_o", inc_d="gdp_d", vce_robust=TRUE, data=grav_small)
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

bvw <- function(y, dist, x, inc_o, inc_d, vce_robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(vce_robust))
  stopifnot(is.character(y), y %in% colnames(data), length(y) == 1)
  stopifnot(is.character(dist), dist %in% colnames(data), length(dist) == 1)
  stopifnot(is.character(x), all(x %in% colnames(data)))
  stopifnot(is.character(inc_d) | inc_d %in% colnames(data) | length(inc_d) == 1)
  stopifnot(is.character(inc_o) | inc_o %in% colnames(data) | length(inc_o) == 1)
  
  # Discarding unusable observations ----------------------------------------
  d <- data %>% 
    filter_at(vars(!!sym(dist)), any_vars(!!sym(dist) > 0)) %>% 
    filter_at(vars(!!sym(dist)), any_vars(is.finite(!!sym(dist)))) %>% 
    
    filter_at(vars(!!sym(y)), any_vars(!!sym(y) > 0)) %>% 
    filter_at(vars(!!sym(y)), any_vars(is.finite(!!sym(y))))
  
  # Transforming data, logging distances ---------------------------------------
  d <- data %>% 
    mutate(
      dist_log = log(!!sym(dist))
    )
  
  # Transforming data, logging flows -------------------------------------------
  d <- d %>% 
    mutate(
      y_log_bvw = log(
        !!sym(y) / (!!sym(inc_o) * !!sym(inc_d))
      )
    )
  
  # GDP weights ----------------------------------------------------------------
  d <- d %>% 
    group_by(!!sym("iso_o")) %>% 
    mutate(inc_world = sum(!!sym(inc_d), na.rm = TRUE)) %>% 
    ungroup()
  
  # same for inc_o or inc_d as we have a squared dataset
  d <- d %>%
    mutate(
      theta_i = !!sym(inc_o) / !!sym("inc_world"),
      theta_j = !!sym(inc_d) / !!sym("inc_world")
    )
  
  # Multilateral resistance (MR) for distance ----------------------------------
  d <- d %>% 
    group_by(!!sym("iso_o"), add = FALSE) %>% 
    mutate(
      mr_dist_1 = sum(!!sym("theta_j") * !!sym("dist_log"), na.rm = TRUE)
    ) %>% 
    
    group_by(!!sym("iso_d"), add = FALSE) %>% 
    mutate(
      mr_dist_2 = sum(!!sym("theta_i") * !!sym("dist_log"), na.rm = TRUE)
    ) %>% 
    
    ungroup() %>% 
    
    mutate(mr_dist_3 = sum(!!sym("theta_i") * !!sym("theta_j") * !!sym("dist_log"), na.rm = TRUE)) %>% 
    
    mutate(dist_log_mr = !!sym("dist_log") - !!sym("mr_dist_1") - !!sym("mr_dist_2") + !!sym("mr_dist_3"))
  
  # Multilateral resistance (MR) for the other independent variables -----------
  d2 <- d %>% 
    select(!!sym("iso_o"), !!sym("iso_d"), !!sym("theta_j"), !!sym("theta_i"), x) %>% 
    gather(!!sym("key"), !!sym("value"), -!!sym("iso_o"), -!!sym("iso_d"), -!!sym("theta_j"), -!!sym("theta_i")) %>% 
    
    mutate(key = paste0(!!sym("key"), "_mr")) %>% 
    
    group_by(!!sym("iso_o"), !!sym("key")) %>% 
    mutate(mr1 = sum(!!sym("theta_j") * !!sym("value"), na.rm = TRUE)) %>% 
    
    group_by(!!sym("iso_d"), !!sym("key")) %>% 
    mutate(mr2 = sum(!!sym("theta_i") * !!sym("value"), na.rm = TRUE)) %>% 
    
    ungroup() %>% 
    mutate(mr3 = sum(!!sym("theta_i") * !!sym("theta_j") * !!sym("value"), na.rm = TRUE)) %>% 
    
    mutate(value = !!sym("value") - !!sym("mr1") - !!sym("mr2") + !!sym("mr3")) %>% 
    
    select(!!!syms(c("iso_o", "iso_d", "key", "value"))) %>% 
    spread(!!sym("key"), !!sym("value"))
  
  # Model ----------------------------------------------------------------------
  dmodel <- left_join(d, d2, by = c("iso_o", "iso_d")) %>% 
    select(!!sym("y_log_bvw"), ends_with("_mr"))
  
  model_bvw <- stats::lm(y_log_bvw ~ ., data = dmodel)
  
  # Return ---------------------------------------------------------------------
  if (vce_robust == TRUE) {
    return_object_1      <- .robustsummary.lm(model_bvw, robust = TRUE)
    return_object_1$call <- as.formula(model_bvw)
    return(return_object_1)
  }
  
  if (vce_robust == FALSE) {
    return_object_1      <- .robustsummary.lm(model_bvw, robust = FALSE)
    return_object_1$call <- as.formula(model_bvw)
    return(return_object_1)
  }
}
