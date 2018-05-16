library(tidyverse)
library(rlang)

load("~/GitHub/gravity/data/gravity_no_zeros.rdata")
data=as.data.frame(gravity_no_zeros)
y="flow"
dependent_variable = y
dist="distw"
distance=dist
x=c("rta")
#x=c("rta", "comcur", "contig")
inc_o="gdp_o"
inc_d="gdp_d"
lab_o = "iso_o"
lab_d = "iso_d"
vce_robust = TRUE

incomes <- c("gdp_o", "gdp_d")
codes <- c("iso_o", "iso_d")
regressors = c("distw", "rta")
code_destination <- lab_d

# summary -----------------------------------------------------------------

.robustsummary.lm <- function(object, robust = FALSE, ...) {
  qr.lm <- function(x, ...) {
    if (is.null(r <- x$qr)) 
      stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
    r
  }
  
  # add extension for robust standard errors
  if (robust == TRUE) { 
    # save variable that are necessary to calcualte robust sd
    X   <- stats::model.matrix(object)
    u2  <- stats::residuals(object)^2
    XDX <- 0
    
    ## One needs to calculate X'DX. But due to the fact that
    ## D is huge (NxN), it is better to do it with a cycle.
    for (i in 1:nrow(X)) {
      XDX <- XDX + u2[i] * X[i,] %*% t(X[i,])
    }
    
    # inverse(X'X)
    XX1 <- solve(t(X) %*% X, tol = 1e-100)
    
    # Sandwich Variance calculation (Bread x meat x Bread)
    varcovar <- XX1 %*% XDX %*% XX1
    
    # adjust degrees of freedom 
    dfc_r <- sqrt(nrow(X)) / sqrt(nrow(X) - ncol(X))
    
    # Standard errors of the coefficient estimates are the
    # square roots of the diagonal elements
    rstdh <- dfc_r * sqrt(diag(varcovar))
  }
  # add extension for clustered standard errors
  
  z   <- object
  p   <- z$rank
  rdf <- z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r   <- sqrt(w) * r
    }
    resvar           <- rss/rdf
    ans              <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans)       <- "summary.lm"
    ans$aliased      <- is.na(stats::coef(object))
    ans$residuals    <- r
    ans$df           <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(
      NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms)) 
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(object, "lm")) 
    warning("calling summary.lm(<fake-lm-object>) ...")
  Qr <- qr.lm(object)
  n  <- NROW(Qr$qr)
  if (is.na(z$df.residual) || n - p != z$df.residual) 
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept")) 
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r   <- sqrt(w) * r
  }
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + stats::var(f)) * 
      1e-30) 
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- 1L:p
  R  <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  
  if (robust == TRUE) se <- rstdh
  
  est  <- z$coefficients[Qr$pivot[p1]]
  tval <- est/se
  ans  <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  pval <- 2 * stats::pt(abs(tval), rdf, lower.tail = FALSE)
  ans$coefficients <- cbind(est, se, tval, pval)
  dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans$aliased <- is.na(stats::coef(object))
  ans$sigma   <- sqrt(resvar)
  ans$df      <- c(p, rdf, NCOL(Qr$qr))
  
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 
      1L
    else 0L
    ans$r.squared     <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - 
                                                       df.int)/rdf)
    ans$fstatistic    <- c(value = (mss/(p - df.int))/resvar, 
                           numdf = p - df.int, dendf = rdf)
    if (robust == TRUE) {
      if (df.int == 0) {
        pos_coef <- 1:length(z$coefficients)
      } else {
        pos_coef <- match(names(z$coefficients)[-match("(Intercept)",
                                                       names(z$coefficients))],
                          names(z$coefficients))
      }
      P_m <- matrix(z$coefficients[pos_coef])
      
      R_m <- diag(1, 
                  length(pos_coef), 
                  length(pos_coef))
      
      ans$fstatistic <- c(value = t(R_m %*% P_m) %*%
                            (solve(varcovar[pos_coef,pos_coef],tol = 1e-100)) %*%
                            (R_m %*% P_m)/(p - df.int), 
                          numdf = p - df.int, dendf = rdf)
      
    }
    
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]
  if (!is.null(z$na.action)) 
    ans$na.action <- z$na.action
  class(ans) <- "summary.lm"
  ans
}
robust_summary <- function(object, robust = FALSE, ...) {
  # Check -------------------------------------------------------------------
  qr_lm <- function(x, ...) {
    if (is.null(r <- x$qr))
      stop("lm object does not have a proper 'qr' component.\n 
           Rank zero or should not have used lm(.., qr=FALSE).")
    r
  }
  
  # Robust standard errors --------------------------------------------------
  if (robust == TRUE) {
    # Save variables that are necessary to calculate robust sd
    X   <- stats::model.matrix(object)
    u2  <- stats::residuals(object)^2
    XDX <- t(X) %*% (X * u2)
    
    # Inverse(X'X)
    XX1 <- solve(t(X) %*% X, tol = 1e-100)
    
    # Sandwich Variance calculation (Bread x meat x Bread)
    varcovar <- XX1 %*% XDX %*% XX1
    
    # Adjust degrees of freedom 
    dfc_r <- sqrt(nrow(X)) / sqrt(nrow(X) - ncol(X))
    
    # Standard errors of the coefficient estimates are the
    # square roots of the diagonal elements
    rstdh <- dfc_r * sqrt(diag(varcovar))
  }
  
  # Clustered standard errors -----------------------------------------------
  z   <- object
  p   <- z$rank
  rdf <- z$df.residual
  
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r   <- sqrt(w) * r
    }
    resvar           <- rss / rdf
    ans              <- z[c("call", "terms", 
                            if (!is.null(z$weights)) "weights")]
    class(ans)       <- "summary.lm"
    ans$aliased      <- is.na(stats::coef(object))
    ans$residuals    <- r
    ans$df           <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(
      NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  
  if (is.null(z$terms)) 
    stop("invalid 'lm' object:  no 'terms' component")
  
  if (!inherits(object, "lm")) 
    warning("calling summary.lm(<fake-lm-object>) ...")
  
  Qr <- qr_lm(object)
  n  <- NROW(Qr$qr)
  
  if (is.na(z$df.residual) || n - p != z$df.residual) 
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept")) 
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r   <- sqrt(w) * r
  }
  
  resvar <- rss / rdf
  
  if (is.finite(resvar) && resvar < (mean(f)^2 + stats::var(f)) * 1e-30)
    warning("essentially perfect fit: summary may be unreliable")
  
  p1 <- 1L:p
  R  <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  
  if (robust == TRUE) se <- rstdh
  
  est  <- z$coefficients[Qr$pivot[p1]]
  tval <- est / se
  ans  <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  pval <- 2 * stats::pt(abs(tval), rdf, lower.tail = FALSE)
  ans$coefficients <- cbind(est, se, tval, pval)
  dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
                                     c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  ans$aliased <- is.na(stats::coef(object))
  ans$sigma   <- sqrt(resvar)
  ans$df      <- c(p, rdf, NCOL(Qr$qr))
  
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept")) 
      1L
    else 0L
    ans$r.squared     <- mss / (mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int) / rdf)
    ans$fstatistic    <- c(value = (mss/(p - df.int)) / resvar, 
                           numdf = p - df.int, dendf = rdf)
    if (robust == TRUE) {
      if (df.int == 0) {
        pos_coef <- 1:length(z$coefficients)
      } else {
        pos_coef <- match(names(z$coefficients)[-match("(Intercept)", 
                                                       names(z$coefficients))],
                          names(z$coefficients))
      }
      P_m <- matrix(z$coefficients[pos_coef])
      
      R_m <- diag(1, 
                  length(pos_coef), 
                  length(pos_coef))
      
      ans$fstatistic <- c(value = t(R_m %*% P_m) %*%
                            (solve(varcovar[pos_coef,pos_coef], tol = 1e-100)) %*%
                            (R_m %*% P_m) / (p - df.int), 
                          numdf = p - df.int, dendf = rdf)
      
    }
  }
  else {
    ans$r.squared <- 0
    ans$adj.r.squared <- 0
  }
  
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,1)]
  
  if (!is.null(z$na.action)) 
    ans$na.action <- z$na.action
  
  # Output ------------------------------------------------------------------
  class(ans) <- "summary.lm"
  ans
}

# tetrads original --------------------------------------------------------

multiway_vcov = T
k = "USA"
ell = "JPN"
Tetrads <- function(){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  
  
  if((multiway_vcov %in% c(TRUE, FALSE)) == FALSE)                        stop("'multiway_vcov' has to be either 'TRUE' or 'FALSE'")
  if(!k%in%data$iso_d)                                                    stop("'k' must be in 'data$iso_d'")
  if(!ell%in%data$iso_o)                                                  stop("'ell' must be in 'data$iso_o'")
  
  # Transforming data, log flows and distances ---------------------------------
  
  d           <- data
  d$dist_log  <- (log(d[dist][,1]))
  d$y_log     <- log(d[y][,1])
  
  # Truncating dataset to only those countries which trade with reference
  # importer and exporter ------------------------------------------------------
  
  d_2a        <- d[d$iso_d == k,] # all iso_o which should stay in dataset
  d_2b        <- d_2a$iso_o
  d_2         <- d[d$iso_o %in% d_2b,] 
  d_3a        <- d_2[d_2$iso_o == ell,] # all iso_d which should stay in dataset
  d_3b        <- d_3a$iso_d
  d_3         <- d_2[d_2$iso_d %in% d_3b,] 
  num.ind.var <- length(x) #independent variables apart from distance
  rm(data); rm(d); rm(d_2a);rm(d_2b); rm(d_2); rm(d_3a); rm(d_3b)
  
  # Taking ratios, ratk --------------------------------------------------------
  
  d_3$lXinratk  <- NA
  d_3$ldistratk <- NA
  
  for(i in unique(d_3$iso_o)){
    
    d_3[d_3$iso_o==i,]$lXinratk <- d_3[d_3$iso_o==i,]$y_log - d_3[d_3$iso_o==i & d_3$iso_d==k,]$y_log
    d_3[d_3$iso_o==i,]$ldistratk <- d_3[d_3$iso_o==i,]$dist_log - d_3[d_3$iso_o==i & d_3$iso_d==k,]$dist_log
  }
  
  # Taking ratios, ratk, for the other independent variables -------------------
  
  ind.var.ratk                  <- list(length=num.ind.var+1)
  ind.var.ratk[[num.ind.var+1]] <- d_3$iso_o
  
  for(j in 1:num.ind.var){
    ind.var.ratk[[j]] <- rep(NA, times=nrow(d_3))
  }
  
  for(j in 1:num.ind.var){
    for(i in unique(d_3$iso_o)){
      ind.var.ratk[[j]][ind.var.ratk[[num.ind.var+1]]==i] <- ((d_3[d_3$iso_o==i,])[x[j]])[,1] - (d_3[d_3$iso_o==i & d_3$iso_d==k,])[x[j]][,1]
    }
  }
  
  for(j in 1:length(x)){
    d_3[x[j]]  <-  ind.var.ratk[[j]]
  }
  
  # Taking the ratio of ratios, rat --------------------------------------------
  
  d_3$lXinrat  <- NA
  d_3$ldistrat <- NA
  
  for(i in unique(d_3$iso_d)){
    
    d_3[d_3$iso_d==i,]$lXinrat  <- d_3[d_3$iso_d==i,]$lXinratk - d_3[d_3$iso_d==i & d_3$iso_o==ell,]$lXinratk
    d_3[d_3$iso_d==i,]$ldistrat <- d_3[d_3$iso_d==i,]$ldistratk - d_3[d_3$iso_d==i & d_3$iso_o==ell,]$ldistratk
  }
  
  d_3$y_log_rat    <- d_3$lXinrat
  d_3$dist_log_rat <- d_3$ldistrat
  
  # Taking the ratio of ratios, rat, for the other independent variables -------
  
  ind.var.rat                  <- list(length=num.ind.var+1)
  ind.var.rat[[num.ind.var+1]] <- d_3$iso_d
  
  for(j in 1:num.ind.var){
    ind.var.rat[[j]] <- rep(NA, times=nrow(d_3))
  }
  
  for(j in 1:num.ind.var){
    for(i in unique(d_3$iso_d)){
      
      ind.var.rat[[j]][ind.var.rat[[num.ind.var+1]]==i] <- ((d_3[d_3$iso_d==i,])[x[j]])[,1] - (d_3[d_3$iso_d==i & d_3$iso_o==ell,])[x[j]][,1]
    }
  }
  
  for(j in 1:length(x)){
    d_3[x[j]]  <-  ind.var.rat[[j]]
  }
  
  # Model ----------------------------------------------------------------------
  
  x_rat <- paste0(x,"_rat")
  
  # new row in dataset for independent _rat variable
  for(j in x){
    l        <- which(x == j)
    rat      <- x_rat[l]
    d_3[rat] <- NA
    d_3[rat] <- d_3[x[l]]
  }
  
  vars  <- paste(c("dist_log_rat", x_rat), collapse=" + ")
  form  <- paste("y_log_rat","~",vars)
  form2 <- stats::as.formula(form)
  
  model_tetrads        <- stats::lm(form2, data = d_3)
  cluster.formula      <- ~ iso_o + iso_d
  model_tetrads_vcov   <- multiwayvcov::cluster.vcov(model=model_tetrads, cluster = cluster.formula)
  model_tetrads_robust <- lmtest::coeftest(x=model_tetrads, vcov=model_tetrads_vcov)
  
  # Return ---------------------------------------------------------------------
  
  if(multiway_vcov == TRUE){
    summary_ted              <- .robustsummary.lm(model_tetrads, robust=TRUE)
    summary_ted$coefficients <- model_tetrads_robust[1:length(rownames(model_tetrads_robust)),]
    return_object            <- summary_ted
    return_object$call       <- form2
    return(return_object)}
  
  if(multiway_vcov == FALSE){
    return_object            <- .robustsummary.lm(model_tetrads, robust=FALSE)
    return_object$call       <- form2
    return(return_object)}
  
}
Tetrads()

# tetrads nueva -----------------------------------------------------------

reference_countries = c("JPN", "USA")
multiway = TRUE
#comente 679 y 685 sino da error pq no encuentra form2
tetrads <- function() {
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
tetrads()
