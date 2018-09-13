context("test-gravity.R")

test_that("get_data connects to the API and returns a valid tibble after valid input", {
  # fit model with example dataset
  fit <- ppml(dependent_variable = "flow", regressors = c("distw", "rta","iso_o","iso_d"),
              robust = TRUE, data = gravity::gravity_zeros)
  
  expect_is(fit, "summary.lm")
  expect_is(fit$coefficients, "matrix")
  expect_output(str(fit), "List of 9")
})
