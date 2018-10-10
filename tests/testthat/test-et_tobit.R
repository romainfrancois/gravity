context("test-gravity.R")

test_that("ET Tobit returns a valid output", {
  # fit model with example dataset
  data("gravity_no_zeros")
  countries_chosen <- names(sort(table(gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
  grav_small <- gravity_no_zeros[gravity_no_zeros$iso_o %in% countries_chosen, ]

  library(dplyr)

  grav_small <- grav_small %>%
    mutate(
      flow = ifelse(flow < 5, 0, flow), # just for testing purposes
      lgdp_o = log(gdp_o),
      lgdp_d = log(gdp_d)
    )

  fit <- et_tobit(
    dependent_variable = "flow", distance = "distw", additional_regressors = c("rta", "lgdp_o", "lgdp_d"),
    data = grav_small
  )

  expect_is(fit, "censReg")
})
