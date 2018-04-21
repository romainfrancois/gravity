## ----read, eval = FALSE--------------------------------------------------
#  library(foreign)
#  
#  col_regfile09 <- read.dta("col_regfile09.dta")

## ----isolate, eval = F---------------------------------------------------
#  data06 <- col_regfile09[col_regfile09$year == 2006,]

## ----choose, eval= F-----------------------------------------------------
#  data06 <- data06[,c(2, 3, 6, 8, 12, 27, 34, 4, 5, 29)]

## ----complete-cases, eval = F--------------------------------------------
#  data06 <- data06[complete.cases(data06) == TRUE,]
#  
#  gravity_zeros <- data06
#  
#  row.names(gravity_zeros) <- 1:length(row.names(gravity_zeros))

## ----scaling, eval = F---------------------------------------------------
#  gravity_zeros$gdp_o <- gravity_zeros$gdp_o / 1000000
#  
#  gravity_zeros$gdp_d <- gravity_zeros$gdp_d / 1000000

