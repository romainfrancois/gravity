context("test-gravity.R")

test_that("BVU returns a valid output", {
  # fit model with example dataset
  data("gravity_no_zeros")
  countries_chosen <- names(sort(table(gravity_no_zeros$iso_o), decreasing = TRUE)[1:10])
  grav_small <- gravity_no_zeros[gravity_no_zeros$iso_o %in% countries_chosen, ]

  fit <- bvu(
    dependent_variable = "flow", distance = "distw", additional_regressors = c("rta", "contig", "comcur"),
    income_origin = "gdp_o", income_destination = "gdp_d", code_origin = "iso_o", code_destination = "iso_d",
    robust = TRUE, data = grav_small
  )

  expect_is(fit, "lm")
})
