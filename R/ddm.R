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
#' @param y name (type: character) of the dependent variable in the dataset 
#' \code{data}, e.g. trade flows. This variable is logged and then used as the 
#' dependent variable in the estimation.
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
#' Unilateral effects drop out due to double demeaning and therefore 
#' cannot be estimated.
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
#' ddm(y="flow", dist="distw", x=c("rta"), 
#' vce_robust=TRUE, data=gravity_no_zeros)
#' 
#' ddm(y="flow", dist="distw", x=c("rta", "comcur", "contig"), 
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
#' ddm(y="flow", dist="distw", x=c("rta"), vce_robust=TRUE, data=grav_small)
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

ddm <- function(y, dist, x, vce_robust=TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(vce_robust))
  stopifnot(is.character(y), y %in% colnames(data), length(y) == 1)
  stopifnot(is.character(dist), dist %in% colnames(data), length(dist) == 1)
  stopifnot(is.character(x), all(x %in% colnames(data)))

  # Transforming data, logging distances ---------------------------------------
  d <- data
  d <- d %>% 
    mutate(
      dist_log_dd = log(!!sym(dist))
    )
  
  # Transforming data, logging flows -------------------------------------------
  d <- d %>% 
    mutate(
      y_log_dd = log(!!sym(y))
    )
  
  # Substracting the means -----------------------------------------------------
  d <- d %>% 
    group_by(!!sym("iso_o"), add = FALSE) %>% 
    mutate(
      y_log_dd = y_log_dd - mean(y_log_dd),
      dist_log_dd = dist_log_dd - mean(dist_log_dd)
    ) %>% 
    
    group_by(!!sym("iso_d"), add = FALSE) %>% 
    mutate(
      y_log_dd = y_log_dd - mean(y_log_dd),
      dist_log_dd = dist_log_dd - mean(dist_log_dd)
    ) %>% 
    
    ungroup() %>% 
    mutate(
      y_log_dd = y_log_dd + mean(y_log_dd),
      dist_log_dd = dist_log_dd + mean(dist_log_dd)
    )
  
  # Substracting the means for the other independent variables -----------------
  d2 <- d %>% 
    select(!!sym("iso_o"), !!sym("iso_d"), x) %>% 
    gather(!!sym("key"), !!sym("value"), -!!sym("iso_o"), -!!sym("iso_d")) %>% 
    
    group_by(!!sym("iso_o"), !!sym("key")) %>% 
    mutate(mean.ind.var.1 = mean(!!sym("value"))) %>% 
    
    group_by(!!sym("iso_d"), !!sym("key"), add = FALSE) %>% 
    mutate(mean.ind.var.2 = mean(!!sym("value"))) %>% 
    
    group_by(!!sym("iso_o"), !!sym("key"), add = FALSE) %>% 
    mutate(dd1 = !!sym("value") - mean.ind.var.1) %>% 
    
    group_by(!!sym("iso_d"), add = FALSE) %>% 
    mutate(dd1 = dd1 - mean.ind.var.2) %>% 
    
    ungroup() %>% 
    mutate(dd1 = dd1 + mean(dd1))
  
  
  for (j in 1:length(x)) {
    ind.var.dd[[j]]     <- d[x[j]][,1]
    mean.ind.var.1[[j]] <- tapply(d[x[j]][,1], d$iso_o, mean)
    mean.ind.var.2[[j]] <- tapply(d[x[j]][,1], d$iso_d, mean)
  }
  
  w   <- letters[1:length(x)]
  
  d_2 <- d
  for (j in 1:length(x)) {
    d_2[w[j]] <- ind.var.dd[[j]]
  }

  d_3 <- d_2
  
  for (j in 1:length(x)) {
    
    for (i in unique(d_2$iso_o)) {
      d_2[d_2$iso_o == i,][w[j]] <- d_2[d_2$iso_o == i,][w[j]] - mean.ind.var.1[[j]][i]
    }
    
    for (i in unique(d_2$iso_d)) {
      d_2[d_2$iso_d == i,][w[j]] <- d_2[d_2$iso_d == i,][w[j]] - mean.ind.var.2[[j]][i]
    }
    
    d_2[w[j]] <- d_2[w[j]] + mean(d_2[x[j]][,1])
    d_3[x[j]] <- d_2[w[j]]
  }
  
  # Model ----------------------------------------------------------------------
  
  x_dd <- paste0(x,"_dd")
  
  # new row in dataset for independent _dd variable
  for (j in x) {
    l       <- which(x == j)
    dd      <- x_dd[l]
    d_3[dd] <- NA
    d_3[dd] <- d_3[x[l]]
  }
  
  vars      <- paste(c("dist_log_dd", x_dd), collapse = " + ")
  form      <- paste("y_log_dd", "~", vars, "- 1")
  form2     <- stats::as.formula(form)
  
  model.ddm <- stats::lm(form2, data = d_3)   

  # Return ---------------------------------------------------------------------
  
  if (vce_robust == TRUE) {
    return.object.1         <- .robustsummary.lm(model.ddm, robust = TRUE)
    return.object.1$call    <- form2
    return(return.object.1)
  }
  
  if (vce_robust == FALSE) {
    return.object.1        <- .robustsummary.lm(model.ddm, robust = FALSE)
    return.object.1$call   <- form2
    return(return.object.1)
  }
}
