context("test-gravity.R")

test_that("get_data connects to the API and returns a valid tibble after valid input", {
  # fit model with example dataset
  library(dplyr)
  
  gravity_no_zeros <- gravity::gravity_no_zeros %>% 
    mutate(
      lgdp_o = log(gdp_o),
      lgdp_d = log(gdp_d)
    )
  
  fit <- gpml(dependent_variable = "flow", regressors = c("distw", "rta", "lgdp_o", "lgdp_d"),
              robust = TRUE, data = gravity_no_zeros)

  expect_is(fit, "summary.lm")
  expect_is(fit$coefficients, "matrix")
  expect_output(str(fit), "List of 9")
})
