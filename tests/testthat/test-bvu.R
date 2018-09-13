context("test-gravity.R")

test_that("get_data connects to the API and returns a valid tibble after valid input", {
  # fit model with example dataset
  fit <- bvu(dependent_variable = "flow", regressors = c("distw", "rta", "contig", "comcur"),
             incomes = c("gdp_o", "gdp_d"), codes = c("iso_o", "iso_d"),
             robust = TRUE, data = gravity::gravity_no_zeros)
  
  expect_is(fit, "summary.lm")
  expect_is(fit$coefficients, "matrix")
  expect_output(str(fit), "List of 11")
})
