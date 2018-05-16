#' @title Tetrads
#' 
#' @description \code{tetrads} estimates gravity models 
#' by taking the ratio of the ratio of flows.
#' 
#' @details \code{tetrads} is an estimation method for gravity models
#' developed by \insertCite{Head2010;textual}{gravity}.
#' 
#' The function \code{tetrads} utilizes the multiplicative form of the
#' gravity equation. After choosing a reference importer \code{K} and 
#' exporter \code{L} one can eliminate importer and exporter fixed effects 
#' by taking the ratio of ratios. 
#' 
#' Only those exporters trading with the 
#' reference importer and importers trading with the reference exporter will 
#' remain for the estimation. Therefore, reference countries should
#' preferably be countries which trade with every other country in the dataset. 
#' 
#' After restricting the data in this way, \code{tetrads} estimates the gravity 
#' equation in its additive form by OLS.
#' 
#' By taking the ratio of ratios, all monadic effects diminish, hence no
#' unilateral variables such as GDP can be included as independent variables.
#' 
#' \code{tetrads} estimation can be used for both, cross-sectional as well as 
#' panel data. Nonetheless, the function is designed to be consistent with the 
#' Stata code for cross-sectional data provided on the website
#' \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook}
#' when choosing robust estimation.
#' 
#' The function \code{tetrads} was therefore tested for cross-sectional data.
#' 
#' If tetrads is used for panel data, the user may has to drop distance as an
#' independent variable as time-invariant effects drop.
#' 
#' For applying \code{tetrads} to panel data see \insertCite{Head2010;textual}{gravity}.
#' 
#' @param dependent_variable name (type: character) of the dependent variable in the dataset 
#' \code{data} (e.g. trade flows).
#' 
#' This variable is logged and then used as the dependent variable in the estimation. 
#' 
#' @param regressors name (type: character) of the regressors to include in the model.
#' 
#' Include the distance variable in the dataset \code{data} containing a measure of 
#' distance between all pairs of bilateral partners and bilateral variables that should 
#' be taken as the independent variables in the estimation. 
#' 
#' The distance is logged automatically when the function is executed.
#' 
#' Unilateral effects drop as the ratio of ratios is taken.
#' 
#' Write this argument as \code{c(distance, contiguity, common curreny, ...)}.
#' 
#' @param codes variable name (type: character) of the code of the country 
#' of origin and destination (e.g. ISO-3 codes from the variables \code{iso_o} and \code{iso_d}) in the 
#' example datasets). 
#' 
#' The variables are grouped by using \code{iso_o} and \code{iso_d} to obtain estimates.
#' 
#' Write this argument as \code{c(code origin, code destination)}.
#' 
#' @param reference_countries reference exporting and importing country, default is set to 
#' \code{c("JPN", "USA")}
#' 
#' Write this argument as \code{c(importing country, exporting country)}.
#' 
#' @param multiway (type: logic) optional; determines whether a function
#' implementing \insertCite{Cameron2011;textual}{gravity} multi-way clustering of 
#' variance-covariance matrices in the package \code{\link[multiwayvcov]} is used
#' for the estimation. 
#' 
#' In case \code{multiway = TRUE}, the \code{\link[multiwayvcov]{cluster.vcov}} \function is used. 
#' 
#' The default value is set to \code{TRUE}.
#'  
#' @param data name of the dataset to be used (type: character).
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
#' When using panel data, a variable for the time may be included in the 
#' dataset. Note that the variable for the time dimension should be of 
#' type factor.
#' 
#' The time variable can be used as a single dependent variable or interaction 
#' term with other variables such as country identifiers by inserting it into 
#' \code{regressors} or as an optional parameter.
#' 
#' The function will remove zero flows and distances.
#' 
#' @param ... additional arguments to be passed to functions used by 
#' \code{tetrads}.
#' 
#' @references 
#' For information on \code{tetrads} see
#' 
#' \insertRef{Cameron2011}{gravity}
#' 
#' \insertRef{Head2010}{gravity}
#' 
#' For more information on gravity models, theoretical foundations and
#' estimation methods in general see 
#' 
#' \insertRef{Anderson1979}{gravity}
#' 
#' \insertRef{Anderson2001}{gravity}
#' 
#' \insertRef{Anderson2010}{gravity}
#' 
#' \insertRef{Head2010}{gravity}
#' 
#' \insertRef{Head2014}{gravity}
#' 
#' \insertRef{Santos2006}{gravity}
#' 
#' and the citations therein.
#' 
#' See \href{https://sites.google.com/site/hiegravity/}{Gravity Equations: Workhorse, Toolkit, and Cookbook} for gravity datasets and Stata code for estimating gravity models.
#' 
#' For estimating gravity equations using panel data see 
#' 
#' \insertRef{Egger2003}{gravity}
#' 
#' \insertRef{Gomez-Herrera2013}{gravity}
#' 
#' and the references therein.
#' 
#' @examples 
#' \dontrun{
#' data(gravity_no_zeros)
#' 
#' tetrads(dependent_variable = "flow", regressors = c("distw", "rta"), 
#' codes = c("iso_o", "iso_d"), reference_countries = c("JPN", "USA"), 
#' multiway = TRUE, data = gravity_no_zeros)
#' 
#' tetrads(dependent_variable = "flow", regressors = c("distw", "rta", "comcur", "contig"), 
#' codes = c("iso_o", "iso_d"), reference_countries = c("JPN", "USA"), 
#' multiway = FALSE, data = gravity_no_zeros)
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
#' tetrads(dependent_variable ="flow", 
#' regressors = c("distw", "rta"), 
#' codes = c("iso_o", "iso_d"), 
#' reference_countries = c(countries_chosen[1], countries_chosen[2]), 
#' multiway = FALSE, 
#' data = grav_small)
#' }
#' 
#' @return
#' The function returns the summary of the estimated gravity model as an 
#' \code{\link[stats]{lm}}-object.
#' 
#' @seealso \code{\link[stats]{lm}}, \code{\link[lmtest]{coeftest}}, 
#' \code{\link[multiwayvcov]{cluster.vcov}}
#' 
#' @export 

tetrads <- function(dependent_variable, regressors, codes, reference_countries = c("JPN", "USA"), 
                    multiway = TRUE, data, ...) {
  # Checks ------------------------------------------------------------------
  stopifnot(is.data.frame(data))
  stopifnot(is.logical(multiway))
  stopifnot(is.character(dependent_variable), dependent_variable %in% colnames(data), length(dependent_variable) == 1)
  stopifnot(is.character(regressors), all(regressors %in% colnames(data)), length(regressors) > 1)
  stopifnot(is.character(reference_countries) | all(reference_countries %in% colnames(data)) | length(reference_countries) == 2)
  
  # Split input vectors -----------------------------------------------------
  code_o <- codes[1]
  code_d <- codes[2]
  
  filter_o <- reference_countries[1]
  filter_d <- reference_countries[2]
  
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
      y_log_tetrads = log(!!sym(dependent_variable))
    )
  
  # Truncating dataset to only those countries which trade with reference
  # importer and exporter ------------------------------------------------------
  d2_filter_d <- d %>% 
    filter_at(vars(!!sym(code_d)), any_vars(!!sym(code_d) == filter_d)) %>% 
    select(!!sym(code_o)) %>% 
    distinct() %>% 
    as_vector()
  
  d2 <- d %>% 
    filter_at(vars(!!sym(code_o)), any_vars(!!sym(code_o) %in% d2_filter_d))
  
  d2_filter_o <- d2 %>% 
    filter_at(vars(!!sym(code_o)), any_vars(!!sym(code_o) == filter_o)) %>% 
    select(!!sym(code_d)) %>% 
    distinct() %>% 
    as_vector()
  
  d2 <- d2 %>% 
    filter_at(vars(!!sym(code_d)), any_vars(!!sym(code_d) %in% d2_filter_o))
  
  # Taking ratios, ratk --------------------------------------------------------
  
  # DESDE ACÁ NO SE COMO HACER UN "NESTED GROUP BY"
  # d_3 <- d2
  
  d2 <- d2 %>% 
    filter(iso_d == filter_d) %>% 
    select(iso_o, y_log_tetrads_d = y_log_tetrads, dist_log_d = dist_log) %>% 
    left_join(d2, ., by = "iso_o") %>% 
    mutate(
      lXinratk = y_log_tetrads - y_log_tetrads_d,
      ldistratk = dist_log - dist_log_d
    ) %>% 
    select(-y_log_tetrads_d, -dist_log_d)
  
  # d_3$lXinratk  <- NA
  # d_3$ldistratk <- NA
  # 
  # for (i in unique(d_3$iso_o)) {
  #   d_3[d_3$iso_o == i,]$lXinratk <- d_3[d_3$iso_o == i,]$y_log_tetrads - d_3[d_3$iso_o == i & d_3$iso_d == filter_d,]$y_log_tetrads
  #   d_3[d_3$iso_o == i,]$ldistratk <- d_3[d_3$iso_o == i,]$dist_log - d_3[d_3$iso_o == i & d_3$iso_d == filter_d,]$dist_log
  # }
  
  # Taking ratios, ratk, for the other independent variables -------------------
  d2 <- d2 %>% 
    select(c("iso_o", "iso_d", regressors)) %>% 
    gather(key, value, -iso_o, -iso_d) %>% 
    left_join(
      d2 %>% 
        filter(iso_d == filter_d) %>% 
        select(c("iso_o", regressors)) %>% 
        gather(key, value, -iso_o)    ,
      by = c("iso_o", "key")
    ) %>% 
    mutate(value = value.x - value.y) %>% 
    select(c("iso_o", "iso_d", "key", "value")) %>% 
    spread(key, value) %>% 
    left_join(d2 %>% select(-one_of(regressors)), by = c("iso_o", "iso_d"))
  
  # num.ind.var <- length(regressors)
  # ind.var.ratk <- list(length = num.ind.var + 1)
  # ind.var.ratk[[num.ind.var + 1]] <- d_3$iso_o
  # 
  # for (j in 1:num.ind.var) {
  #   ind.var.ratk[[j]] <- rep(NA, times = nrow(d_3))
  # }
  # 
  # for (j in 1:num.ind.var) {
  #   for (i in unique(d_3$iso_o)) {
  #     ind.var.ratk[[j]][ind.var.ratk[[num.ind.var + 1]] == i] <- ((d_3[d_3$iso_o == i,])[regressors[j]])[,1][[1]] - (d_3[d_3$iso_o == i & d_3$iso_d == filter_d,])[regressors[j]][,1][[1]]
  #   }
  # }
  # 
  # for (j in 1:length(regressors)) {
  #   d_3[regressors[j]]  <-  ind.var.ratk[[j]]
  # }
  
  # Taking the ratio of ratios, rat --------------------------------------------
  d2 <- d2 %>% 
    filter(iso_o == filter_o) %>% 
    select(iso_d, lXinratk_o = lXinratk, ldistratk_o = ldistratk) %>% 
    left_join(d2, ., by = "iso_d") %>% 
    mutate(
      lXinrat = lXinratk - lXinratk_o,
      ldistrat = ldistratk - ldistratk_o
    ) %>% 
    select(-lXinratk_o, -ldistratk_o) %>% 
    mutate(
      y_log_rat = lXinrat,
      dist_log_rat = ldistrat
    )
  
  # d_3$lXinrat  <- NA
  # d_3$ldistrat <- NA
  # 
  # for (i in unique(d_3$iso_d)) {
  #   d_3[d_3$iso_d == i,]$lXinrat  <- d_3[d_3$iso_d == i,]$lXinratk - d_3[d_3$iso_d == i & d_3$iso_o == filter_o,]$lXinratk
  #   d_3[d_3$iso_d == i,]$ldistrat <- d_3[d_3$iso_d == i,]$ldistratk - d_3[d_3$iso_d == i & d_3$iso_o == filter_o,]$ldistratk
  # }
  # 
  # d_3$y_log_rat    <- d_3$lXinrat
  # d_3$dist_log_rat <- d_3$ldistrat
  
  # Taking the ratio of ratios, rat, for the other independent variables -------
  # ind.var.rat <- list(length = num.ind.var + 1)
  # ind.var.rat[[num.ind.var + 1]] <- d_3$iso_d
  # 
  # for (j in 1:num.ind.var) {
  #   ind.var.rat[[j]] <- rep(NA, times = nrow(d_3))
  # }
  # 
  # for (j in 1:num.ind.var) {
  #   for (i in unique(d_3$iso_d)) {
  #     ind.var.rat[[j]][ind.var.rat[[num.ind.var + 1]] == i] <- ((d_3[d_3$iso_d == i,])[regressors[j]])[,1][[1]] - (d_3[d_3$iso_d == i & d_3$iso_o == filter_o,])[regressors[j]][,1][[1]]
  #   }
  # }
  # 
  # for (j in 1:length(regressors)) {
  #   d_3[regressors[j]]  <-  ind.var.rat[[j]]
  # }
  
  d2 <- d2 %>% 
    select(c("iso_o", "iso_d", additional_regressors)) %>% 
    gather(key, value, -iso_o, -iso_d) %>% 
    left_join(
      d2 %>% 
        filter(iso_o == filter_o) %>% 
        select(c("iso_d", additional_regressors)) %>% 
        gather(key, value, -iso_d)    ,
      by = c("iso_d", "key")
    ) %>% 
    mutate(
      key = paste0(key, "_rat"),
      value = value.x - value.y
    ) %>% 
    select(c("iso_o", "iso_d", "key", "value")) %>% 
    spread(key, value) %>% 
    left_join(d2 %>% select(-one_of(regressors)), by = c("iso_o", "iso_d")) %>% 
    select(ends_with("_rat"), codes)
  
  # Model ----------------------------------------------------------------------
  # x_rat <- paste0(regressors,"_rat")
  # 
  # # new row in dataset for independent _rat variable
  # for (j in regressors) {
  #   l        <- which(regressors == j)
  #   rat      <- x_rat[l]
  #   d_3[rat] <- NA
  #   d_3[rat] <- d_3[regressors[l]]
  # }
  
  additional_regressors <- paste0(additional_regressors, "_rat")
  form <- stats::as.formula(paste("y_log_rat", "~ dist_log_rat +", paste(additional_regressors, collapse = " + ")))
  model_tetrads <- stats::lm(form, data = d2)
  model_tetrads_vcov <- multiwayvcov::cluster.vcov(model = model_tetrads, cluster = ~ iso_o + iso_d)
  model_tetrads_robust <- lmtest::coeftest(x = model_tetrads, vcov = model_tetrads_vcov)
  
  # vars  <- paste(c("dist_log_rat", x_rat), collapse = " + ")
  # form  <- paste("y_log_rat","~", vars)
  # form2 <- stats::as.formula(form)
  # 
  # model_tetrads        <- stats::lm(form2, data = d_3)
  # cluster.formula      <- ~ iso_o + iso_d
  # model_tetrads_vcov   <- multiwayvcov::cluster.vcov(model = model_tetrads, cluster = cluster.formula)
  # model_tetrads_robust <- lmtest::coeftest(x = model_tetrads, vcov = model_tetrads_vcov)
  
  # Return ---------------------------------------------------------------------
  if (multiway == TRUE) {
    summary_ted              <- robust_summary(model_tetrads, robust = TRUE)
    summary_ted$coefficients <- model_tetrads_robust[1:length(rownames(model_tetrads_robust)),]
    return_object            <- summary_ted
    #return_object$call       <- form2
    return(return_object)
  }
  
  if (multiway == FALSE) {
    return_object            <- robust_summary(model_tetrads, robust = FALSE)
    #return_object$call       <- form2
    return(return_object)}
}
