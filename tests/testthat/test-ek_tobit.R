context("test-gravity.R")

test_that("get_data connects to the API and returns a valid tibble after valid input", {
  # fit model with example dataset
  library(dplyr)
  
  gravity_zeros <- gravity::gravity_zeros %>%
    mutate(
      lgdp_o = log(gdp_o),
      lgdp_d = log(gdp_d)
    )
  
  fit <- ek_tobit(dependent_variable = "flow", regressors = c("distw", "rta","lgdp_o","lgdp_d"),
                  code_destination = "iso_d",
                  robust = TRUE, data = gravity_zeros)
  
  expect_is(fit, "summary.survreg")
  expect_is(fit$coefficients, "numeric")
  expect_output(str(fit), "List of 14")
})
