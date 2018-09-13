context("test-gravity.R")

test_that("get_data connects to the API and returns a valid tibble after valid input", {
  # fit model with example dataset
  library(dplyr)
  
  gravity_zeros <- gravity_zeros %>%
    mutate(
      lgdp_o = log(gdp_o),
      lgdp_d = log(gdp_d)
    )
  
  fit <- tobit(dependent_variable = "flow", regressors = c("distw", "rta", "lgdp_o", "lgdp_d"),
               added_constant = 1, data = gravity_zeros)
  
  expect_is(fit, "summary.maxLik")
  expect_is(fit$estimate, "matrix")
  expect_output(str(fit), "List of 11")
})
