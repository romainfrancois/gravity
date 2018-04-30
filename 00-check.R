library(tidyverse)
library(rlang)

data=as.data.frame(gravity_no_zeros)
y="flow"
dist="distw"
x=c("rta")
inc_o="gdp_o"
inc_d="gdp_d"
vce_robust = TRUE


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

# bvumod ---------------------------------------------------------------------

bvumod <- function() {
  # Transforming data, logging distances ---------------------------------------
  d <- data
  d <- d %>% 
    mutate(
      dist_log = log(!!sym(dist)),
      count = row_number()
    )
  
  # Transforming data, logging flows -------------------------------------------
  d <- d %>% 
    mutate(
      y_inc = !!sym(y) / (!!sym(inc_o) * !!sym(inc_d)),
      y_inc_log = log(!!sym("y_inc"))
    )
  
  # Multilateral Resistance (MR) for distance ----------------------------------
  d <- d %>% 
    group_by(!!sym("iso_o")) %>% 
    mutate(mean.dist_log.1 = mean(!!sym("dist_log"))) %>% 
    group_by(!!sym("iso_d"), add = FALSE) %>% 
    mutate(mean.dist_log.2 = mean(!!sym("dist_log"))) %>% 
    ungroup() %>% 
    mutate(
      mean.dist_log.3 = mean(!!sym("dist_log")),
      dist_log_mr = !!sym("dist_log") - 
        (!!sym("mean.dist_log.1") + !!sym("mean.dist_log.2") - !!sym("mean.dist_log.3"))
    )
  
  # Multilateral Resistance (MR) for the other independent variables -----------
  d2 <- d %>% 
    select(!!sym("iso_o"), !!sym("iso_d"), x) %>% 
    gather(!!sym("key"), !!sym("value"), -!!sym("iso_o"), -!!sym("iso_d")) %>% 
    
    group_by(!!sym("iso_o"), !!sym("key")) %>% 
    mutate(mean.dist_log.1 = mean(!!sym("value"))) %>% 
    
    group_by(!!sym("iso_d"), !!sym("key")) %>% 
    mutate(mean.dist_log.2 = mean(!!sym("value"))) %>% 
    
    group_by(!!sym("key")) %>% 
    mutate(
      mean.dist_log.3 = mean(!!sym("value")),
      dist_log_mr = !!sym("value") - (!!sym("mean.dist_log.1") + !!sym("mean.dist_log.2") - !!sym("mean.dist_log.3"))
    ) %>% 
    
    ungroup() %>% 
    mutate(key = paste0(!!sym("key"), "_mr")) %>% 
    
    select(!!!syms(c("iso_o", "iso_d", "key", "dist_log_mr"))) %>% 
    spread(!!sym("key"), !!sym("dist_log_mr"))
  
  
  # Model ----------------------------------------------------------------------
  dmodel <- left_join(d, d2, by = c("iso_o", "iso_d")) %>% 
    select(!!sym("y_inc_log"), ends_with("_mr"))
  
  model.bvu <- stats::lm(y_inc_log ~ ., data = dmodel)
  
  # Return ---------------------------------------------------------------------
  if (vce_robust == TRUE) {
    return.object.1      <- .robustsummary.lm(model.bvu, robust = TRUE)
    return.object.1$call <- as.formula(model.bvu)
    return(return.object.1)}
  
  if (vce_robust == FALSE) {
    return.object.1      <- .robustsummary.lm(model.bvu, robust = FALSE)
    return.object.1$call <- as.formula(model.bvu)
    return(return.object.1)}
}
bvumod()

# bvuorig -----------------------------------------------------------------

bvuorig <- function(){
  
  if(!is.data.frame(data))                                                stop("'data' must be a 'data.frame'")
  if((vce_robust %in% c(TRUE, FALSE)) == FALSE)                           stop("'vce_robust' has to be either 'TRUE' or 'FALSE'")
  if(!is.character(y)     | !y%in%colnames(data)     | length(y)!=1)      stop("'y' must be a character of length 1 and a colname of 'data'")
  if(!is.character(dist)  | !dist%in%colnames(data)  | length(dist)!=1)   stop("'dist' must be a character of length 1 and a colname of 'data'")
  if(!is.character(x)     | !all(x%in%colnames(data)))                    stop("'x' must be a character vector and all x's have to be colnames of 'data'")  
  
  if(!is.character(inc_d) | !inc_d%in%colnames(data) | length(inc_d)!=1)  stop("'inc_d' must be a character of length 1 and a colname of 'data'")
  if(!is.character(inc_o) | !inc_o%in%colnames(data) | length(inc_o)!=1)  stop("'inc_o' must be a character of length 1 and a colname of 'data'")
  
  # Transforming data, logging distances ---------------------------------------
  
  d               <- data
  d$dist_log      <- (log(d[dist][,1]))
  d$count         <- 1:length(d$iso_o)
  
  # Transforming data, logging flows -------------------------------------------
  
  d$y_inc         <- d[y][,1] / (d[inc_o][,1] * d[inc_d][,1])
  d$y_inc_log     <- log(d$y_inc)
  
  # Multilateral Resistance (MR) for distance ----------------------------------
  
  mean.dist_log.1 <- vector(length=length(unique(d$iso_o)))
  mean.dist_log.2 <- vector(length=length(unique(d$iso_o)))  
  
  for(i in unique(d$iso_o)){
    mean.dist_log.1[i] <- mean(d$dist_log[d$iso_o == i])
  }
  
  for(i in unique(d$iso_o)){
    mean.dist_log.2[i] <- mean(d$dist_log[d$iso_d == i])
  }
  
  mean.dist_log.3 <- mean(d$dist_log)
  d$dist_log_mr <- d$dist_log - (mean.dist_log.1[d$iso_o] + 
                                   mean.dist_log.2[d$iso_d] - 
                                   mean.dist_log.3)
  
  # Multilateral Resistance (MR) for the other independent variables -----------
  
  num.ind.var    <- length(x) #independent variables apart from distance
  
  mean.ind.var.1 <- list(length=num.ind.var)
  mean.ind.var.2 <- list(length=num.ind.var)
  mean.ind.var.3 <- list(length=num.ind.var)
  
  for(j in 1:num.ind.var){
    mean.ind.var.1[[j]]        <- rep(NA, times=length(unique(d$iso_o)))
    mean.ind.var.2[[j]]        <- rep(NA, times=length(unique(d$iso_o)))
    names(mean.ind.var.1[[j]]) <- x[j]
    names(mean.ind.var.2[[j]]) <- x[j]
  }
  
  ind.var <- list(length=num.ind.var)
  
  for(j in 1:num.ind.var){
    for(i in unique(d$iso_o)){
      ind.var[[j]] <- ((d[d$iso_o == i,])[x[j]])[,1]
      names(ind.var[[j]]) <- d$iso_d[d$iso_o == i]
      mean.ind.var.1[[j]][i] <- mean(ind.var[[j]])
    }
  }
  
  for(j in 1:num.ind.var){
    for(i in unique(d$iso_o)){
      ind.var[[j]] <- ((d[d$iso_d == i,])[x[j]])[,1]
      names(ind.var[[j]]) <- d$iso_o[d$iso_d == i]
      mean.ind.var.2[[j]][i] <- mean(ind.var[[j]])
    }
  }
  
  for(j in 1:num.ind.var){
    mean.ind.var.3[[j]] <- mean(d[x[[j]]][,1])
  }
  
  # MR 
  d_2 <- d
  for(k in x){
    l <- which(x == k)
    d_2[k] <- d[k][,1] - (mean.ind.var.1[[l]][d$iso_o] + 
                            mean.ind.var.2[[l]][d$iso_d] - 
                            mean.ind.var.3[[l]])
  }
  
  # Model ----------------------------------------------------------------------
  
  x_mr <- paste0(x,"_mr")
  
  # new row in dataset for independent _mr variable
  for(j in x){
    l       <- which(x == j)
    mr      <- x_mr[l]
    d_2[mr] <- NA
    d_2[mr] <- d_2[x[l]]
  }
  
  vars      <- paste(c("dist_log_mr", x_mr), collapse=" + ")
  form      <- paste("y_inc_log","~",vars)
  form2     <- stats::as.formula(form)
  model.BVU <- stats::lm(form2, data = d_2)
  
  # Return ---------------------------------------------------------------------
  
  if(vce_robust == TRUE){
    return.object.1      <- .robustsummary.lm(model.BVU, robust=TRUE)
    return.object.1$call <- form2
    return(return.object.1)}
  
  if(vce_robust == FALSE){
    return.object.1      <- .robustsummary.lm(model.BVU, robust=FALSE)
    return.object.1$call <- form2
    return(return.object.1)}
}
bvuorig()

# bvwmod ---------------------------------------------------------------------

bvwmod <- function() {
  # Transforming data, logging distances ---------------------------------------
  d <- data
  d <- d %>% 
    mutate(
      dist_log = log(!!sym(dist)),
      count = row_number()
    )
  
  # Transforming data, logging flows -------------------------------------------
  d <- d %>% 
    mutate(
      y_inc = !!sym(y) / (!!sym(inc_o) * !!sym(inc_d)),
      y_inc_log = log(!!sym("y_inc"))
    )
  
  # GDP weights ----------------------------------------------------------------
  d <- d %>% 
    group_by(!!sym("iso_o")) %>% 
    mutate(inc_world = sum(!!sym(inc_d))) %>% 
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
      mr.dist.1 = sum(!!sym("theta_j") * !!sym("dist_log"))
    ) %>% 
    
    group_by(!!sym("iso_d"), add = FALSE) %>% 
    mutate(
      mr.dist.2 = sum(!!sym("theta_i") * !!sym("dist_log"))
    ) %>% 
    
    ungroup() %>% 
    
    mutate(mr.dist.3 = sum(!!sym("theta_i") * !!sym("theta_j") * !!sym("dist_log"))) %>% 
    
    mutate(dist_log_mr = !!sym("dist_log") - !!sym("mr.dist.1") - !!sym("mr.dist.2") + !!sym("mr.dist.3"))
  
  # Multilateral resistance (MR) for the other independent variables -----------
  d2 <- d %>% 
    select(!!sym("iso_o"), !!sym("iso_d"), !!sym("theta_j"), !!sym("theta_i"), x) %>% 
    gather(!!sym("key"), !!sym("value"), -!!sym("iso_o"), -!!sym("iso_d"), -!!sym("theta_j"), -!!sym("theta_i")) %>% 
    
    group_by(!!sym("iso_o"), !!sym("key")) %>% 
    mutate(mr.dist.1 = sum(!!sym("theta_j") * !!sym("value"))) %>% 
    
    group_by(!!sym("iso_d"), !!sym("key")) %>% 
    mutate(mr.dist.2 = sum(!!sym("theta_i") * !!sym("value"))) %>% 
    
    group_by(!!sym("key")) %>% 
    mutate(
      mr.dist.3 = !!sym("theta_i") * !!sym("theta_j") * !!sym("value"),
      dist_log_mr = !!sym("value") - !!sym("mr.dist.1") - !!sym("mr.dist.2") + !!sym("mr.dist.3")
    ) %>% 
    
    ungroup() %>% 
    mutate(key = paste0(!!sym("key"), "_mr")) %>% 
    
    select(!!!syms(c("iso_o", "iso_d", "key", "dist_log_mr"))) %>% 
    spread(!!sym("key"), !!sym("dist_log_mr"))
  
  # Model ----------------------------------------------------------------------
  dmodel <- left_join(d, d2, by = c("iso_o", "iso_d")) %>% 
    select(!!sym("y_inc_log"), ends_with("_mr"))
  
  model.bvu <- stats::lm(y_inc_log ~ ., data = dmodel)
  
  # Return ---------------------------------------------------------------------
  if (vce_robust == TRUE) {
    return.object.1      <- .robustsummary.lm(model.bvu, robust = TRUE)
    return.object.1$call <- as.formula(model.bvu)
    return(return.object.1)}
  
  if (vce_robust == FALSE) {
    return.object.1      <- .robustsummary.lm(model.bvu, robust = FALSE)
    return.object.1$call <- as.formula(model.bvu)
    return(return.object.1)}
}
bvwmod()

# bvworig ---------------------------------------------------------------------

bvworig <- function() {
  # Transforming data, logging distances ---------------------------------------
  
  d           <- data
  d$dist_log  <- log(d[dist][,1])
  d$count     <- 1:length(d$iso_o)
  
  # Transforming data, logging flows -------------------------------------------
  
  d$y_inc     <- d[y][,1] / (d[inc_o][,1] * d[inc_d][,1])
  d$y_inc_log <- log(d$y_inc)
  
  # GDP weights ----------------------------------------------------------------
  
  inc_world   <- tapply(d[inc_d][,1], d$iso_o, sum)
  d$inc_world <- as.numeric(inc_world[d$iso_o])
  # same for inc_o or inc_d as we have a squared dataset
  
  d$theta_i   <- d[inc_o][,1] / d$inc_world
  d$theta_j   <- d[inc_d][,1] / d$inc_world
  
  # Multilateral resistance (MR) for distance ----------------------------------
  
  d$mr.dist.1 <- NA
  d$mr.dist.2 <- NA
  
  for(i in names(inc_world)){
    
    d[d$iso_o == i,]$mr.dist.1 <- sum(d[d$iso_o == i,]$theta_j * d[d$iso_o == i,]$dist_log)
    d[d$iso_d == i,]$mr.dist.2 <- sum(d[d$iso_d == i,]$theta_i * d[d$iso_d == i,]$dist_log)
  }
  
  d$mr.dist.3   <- sum(d$theta_i * d$theta_j * d$dist_log)
  d$dist_log_mr <- d$dist_log - d$mr.dist.1 - d$mr.dist.2 + d$mr.dist.3
  
  # Multilateral resistance (MR) for the other independent variables -----------
  
  num.ind.var   <- length(x) #independent variables apart from distance
  d_2           <- d
  
  for(j in 1:num.ind.var){
    
    mr.1 <- noquote(paste(c(noquote(x[j]),noquote(".mr1")),collapse="")) 
    mr.2 <- noquote(paste(c(noquote(x[j]),noquote(".mr2")),collapse=""))  
    mr.3 <- noquote(paste(c(noquote(x[j]),noquote(".mr3")),collapse="")) 
    mr   <- noquote(paste(c(noquote(x[j]),noquote(".mr")),collapse="")) 
    
    d_2[mr.1] <- NA
    d_2[mr.2] <- NA
    d_2[mr.3] <- NA
    d_2[mr]   <- NA
    
    for(i in names(inc_world)){
      
      d_2[d_2$iso_o == i,][mr.1] <- sum(d_2[d_2$iso_o == i,]$theta_j * d_2[d_2$iso_o == i,][x[j]])
      d_2[d_2$iso_d == i,][mr.2] <- sum(d_2[d_2$iso_d == i,]$theta_i * d_2[d_2$iso_d == i,][x[j]])
    }
    
    d_2[mr.3] <- sum(d_2$theta_i * d_2$theta_j * d_2[x[j]])
    d_2[x[j]] <- d_2[x[j]] - d_2[mr.1] - d_2[mr.2] + d_2[mr.3]
    
  }
  
  # Model ----------------------------------------------------------------------
  
  x_mr <- paste0(x,"_mr")
  
  # new row in dataset for independent _mr variable
  for(j in x){
    l       <- which(x == j)
    mr      <- x_mr[l]
    d_2[mr] <- NA
    d_2[mr] <- d_2[x[l]]
  }
  
  vars      <- paste(c("dist_log_mr", x_mr), collapse=" + ")
  form      <- paste("y_inc_log","~",vars)
  form2     <- stats::as.formula(form)
  model.BVW <- stats::lm(form2, data = d_2)
  
  # Return ---------------------------------------------------------------------
  
  if(vce_robust == TRUE){
    return.object.1      <- .robustsummary.lm(model.BVW, robust=TRUE)
    return.object.1$call <- form2
    return(return.object.1)}
  
  if(vce_robust == FALSE){
    return.object.1      <- .robustsummary.lm(model.BVW, robust=FALSE)
    return.object.1$call <- form2
    return(return.object.1)}
}
bvworig()
