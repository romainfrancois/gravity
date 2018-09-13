context("test-gravity.R")

test_that("get_data connects to the API and returns a valid tibble after valid input", {
  # fit model with example dataset
  fit <- sils(dependent_variable = "flow", regressors = c("distw", "rta"),
              incomes = c("gdp_o", "gdp_d"),
              maxloop = 100, dec_places = 4, robust = TRUE, verbose = FALSE,
              data = gravity::gravity_no_zeros)
  
  expect_is(fit, "summary.lm")
  expect_is(fit$coefficients, "matrix")
  expect_output(str(fit), "List of 11")
})
