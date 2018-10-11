context("test-gravity.R")

test_that("Tetrads returns a valid output", {
  # fit model with example dataset
  data("gravity_no_zeros")
  countries_chosen <- names(sort(table(gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
  grav_small <- gravity_no_zeros[gravity_no_zeros$iso_o %in% countries_chosen, ]
  
  fit <- tetrads(
    dependent_variable ="flow",
    distance = "distw",
    additional_regressors = "rta",
    code_origin = "iso_o",
    code_destination = "iso_d",
    filter_origin = countries_chosen[1],
    filter_destination = countries_chosen[2],
    multiway = FALSE,
    data = grav_small
  )
  
  expect_is(fit, "lm")
})
