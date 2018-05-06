#' @title Non-linear Least Squares (NLS)
#' 
#' @description \code{nls} estimates gravity models in their 
#' multiplicative form via Nonlinear Least Squares.
#' 
#' @details \code{nls} is an estimation method for gravity models
#' belonging to generalized linear models.
#' It is estimated via \code{\link[stats]{glm}} using the gaussian distribution and a log-link.
#' As the method may not lead to convergence when poor
#' starting values are used, the linear predictions, fitted values,
#' and estimated coefficients resulting from a 
#' \code{\link[gravity]{ppml}} estimation are used for the arguments 
#' \code{etastart}, \code{mustart}, and \code{start}.
#' To execute the function a square gravity dataset with all pairs of 
#' countries, ISO-codes for the country of origin and destination, a measure of 
#' distance between the bilateral partners as well as all 
#' information that should be considered as dependent an independent 
#' variables is needed. 
#' Missing bilateral flows as well as incomplete rows should be 
#' excluded from the dataset.  
#' Zero trade flows are allowed.
#' For similar functions, utilizing the multiplicative form via the log-link, 
#' but different distributions, see \code{\link[gravity]{ppml}}, \code{\link[gravity]{gpml}}, and \code{\link[gravity]{nbpml}}.
#' 
#' \code{nls} estimation can be used for both, cross-sectional as well as 
#' panel data. 
#' It is up to the user to ensure that the functions can be applied 
#' to panel data. 
#' Depending on the panel dataset and the variables - 
#' specifically the type of fixed effects - 
#' included in the model, it may easily occur that the model is not computable. 
#' Also, note that by including bilateral fixed effects such as country-pair 
#' effects, the coefficients of time-invariant observables such as distance 
#' can no longer be estimated. 
#' Depending on the specific model, the code of the 
#' respective function may has to be changed in order to exclude the distance 
#' variable from the estimation. 
#' At the very least, the user should take special 
#' care with respect to the meaning of the estimated coefficients and variances 
#' as well as the decision about which effects to include in the estimation. 
#' When using panel data, the parameter and variance estimation of the models 
#' may have to be changed accordingly.
#' For a comprehensive overview of gravity models for panel data 
#' see Egger and Pfaffermayr (2003), Gomez-Herrera (2013) and Head, Mayer and 
#' Ries (2010) as well as the references therein. 
#' 
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. 
#' 
#' @param dist name (type: character) of the distance variable in the dataset 
#' \code{data} containing a measure of distance between all pairs of bilateral
#' partners. It is logged automatically when the function is executed. 
#' 
#' @param x vector of names (type: character) of those bilateral variables in 
#' the dataset \code{data} that should be taken as the independent variables 
#' in the estimation. If an independent variable is a dummy variable,
#' it should be of type numeric (0/1) in the dataset. If an independent variable 
#' is defined as a ratio, it should be logged. 
#' Unilateral variables such as country dummies or incomes can be added. 
#' If unilateral metric variables such as GDPs should be used as independent 
#' variables, those variables have to be logged first and the 
#' logged variable can be used in \code{x}.
#' Interaction terms can be added.
#' 
#' @param vce_robust robust (type: logic) determines whether a robust 
#' variance-covariance matrix should be used. The default is set to \code{TRUE}. 
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
#' When using panel data, a variable for the time may be included in the 
#' dataset. Note that the variable for the time dimension should be of 
#' type: factor. See the references for more information on panel data.
#' 
#' @param ... additional arguments to be passed to functions used by 
#' \code{nls}.
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
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' 
#' For estimating gravity equations using panel data see 
#' 
#' Egger, P., & Pfaffermayr, M. (2003) <DOI:10.1007/s001810200146>
#' 
#' Gomez-Herrera, E. (2013) <DOI:10.1007/s00181-012-0576-2>
#' 
#' and the references therein.
#' 
#' @examples 
#' \dontrun{
#' # Example for data with zero trade flows
#' data(gravity_zeros)
#' 
#' gravity_zeros$lgdp_o <- log(gravity_zeros$gdp_o)
#' gravity_zeros$lgdp_d <- log(gravity_zeros$gdp_d)
#' 
#' nls(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), 
#' vce_robust=TRUE, data=gravity_zeros)
#' 
#' # Example for data without zero trade flows
#' data(gravity_no_zeros)
#' 
#' gravity_no_zeros$lgdp_o <- log(gravity_no_zeros$gdp_o)
#' gravity_no_zeros$lgdp_d <- log(gravity_no_zeros$gdp_d)
#' 
#' nls(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), 
#' vce_robust=TRUE, data=gravity_no_zeros)
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
#' nls(y="flow", dist="distw", x=c("rta","lgdp_o","lgdp_d"), vce_robust=TRUE, data=grav_small_zeros)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model similar to a
#' \code{\link[stats]{glm}}-object.
#' 
#' @seealso \code{\link[stats]{glm}}, \code{\link[lmtest]{coeftest}}, 
#' \code{\link[sandwich]{vcovHC}}
#' 
#' @export

nls <- function(dependent_variable, regressors, vce_robust = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(vce_robust))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  
  # Split input vectors -----------------------------------------------------
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
    ) %>% 
    
    select(
      !!sym("dependent_variable"), !!sym("dist_log"), !!sym("additional_regressors")
    )
  
  # Model ----------------------------------------------------------------------
  form <- stats::as.formula(
    sprintf("%s ~ %s + %s", 
            dependent_variable,
            "dist_log",
            paste(additional_regressors, sep = " + "))
  )
  
  # For nls the starting values are retrieved from the resuts of PPML
  model_PPML       <- stats::glm(form, data = d, family = stats::quasipoisson(link = "log"))
  model_PPML_eta   <- model_PPML$linear.predictors
  model_PPML_mu    <- model_PPML$fitted.values
  model_PPML_start <- model_PPML$coefficients
  
  model_nls <- glm(form, data = d, family = stats::gaussian(link = "log"), 
                   control = list(maxit = 200, trace = FALSE),
                   etastart = model_PPML_eta, # linear predictors
                   mustart = model_PPML_mu, # fitted values
                   start = model_PPML_start) # estimated coefficients
  
  model_nls_robust <- lmtest::coeftest(model_nls, vcov = sandwich::vcovHC(model_nls, "HC1"))
  
  # Return --------------------------------------------------------------------- 
  
  if (vce_robust == TRUE) {
    summary_nls                 <- robust_summary(model_nls, robust = TRUE)
    summary_nls$coefficients    <- model_nls_robust[1:length(rownames(model_nls_robust)),]
    return_object               <- summary_nls
    return_object$call          <- form
    return_object$r.squared     <- NULL 
    return_object$adj.r.squared <- NULL
    return_object$fstatistic    <- NULL
    return(return_object)
  }
  
  if (vce_robust == FALSE) {
    return_object      <- summary(model_nls)
    return_object$call <- form
    return(return_object)
  }
}
