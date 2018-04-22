## ----read, eval = F------------------------------------------------------
#  url <- "http://econ.sciences-po.fr/sites/default/files/file/tmayer/data/col_regfile09.zip"
#  zip <- "col_regfile09.zip"
#  
#  if (!file.exists(zip)) { try(download.file(url, zip)) }
#  try(system("7z e col_regfile09.zip"))
#  
#  library(haven)
#  col_regfile09 <- read_dta("col_regfile09.dta")

## ----isolate, eval = F---------------------------------------------------
#  data06 <- col_regfile09[col_regfile09$year == 2006,]

## ----choose, eval= F-----------------------------------------------------
#  library(dplyr)
#  data06 <- data06 %>%
#    select(iso_o, iso_d, distw, gdp_o, gdp_d, rta, flow, contig, comlang_off, comcur)

## ----complete-cases, eval = F--------------------------------------------
#  gravity_zeros <- data06 %>%
#    filter(complete.cases(.))

## ----scaling, eval = F---------------------------------------------------
#  gravity_zeros <- gravity_zeros %>%
#    mutate(
#      gdp_o = gdp_o / 1000000,
#      gdp_d = gdp_d / 1000000,
#    )

## ----no-zeros, eval = F--------------------------------------------------
#  gravity_no_zeros <- gravity_zeros %>%
#    filter(flow != 0)

## ----export, eval = F----------------------------------------------------
#  save(gravity_zeros, file = "gravity_zeros.rdata", compress = "xz")
#  save(gravity_no_zeros, file = "gravity_no_zeros.rdata", compress = "xz")

