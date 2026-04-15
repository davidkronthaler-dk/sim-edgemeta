#-------------------------------------------------------------------------------
# Unit tests for utility functions
#-------------------------------------------------------------------------------

# Error-handling function
test_that("trySim handles error-free and errored expressions", {
  expect_equal(trySim({
    1 + 1
  }), 2)
  expect_equal(trySim({
    stop("fail")
  }), NA_real_)
  expect_equal(trySim({
    stop("fail")
  }, nx = 3), rep(NA_real_, 3))
})

# Coverage of prediction intervals
test_that("coverpi returns expected output", {
  expect_equal(coverpi(c(0, 1), 0.5), 1)
  expect_equal(coverpi(c(0, 1), 2), 0)
  expect_equal(coverpi(c(1, 2), c(1.5, NA)), 1)
  expect_true(is.na(coverpi(c(0, NA), 0.5)))
  expect_true(is.na(coverpi(c(0), 0.5)))
  expect_true(is.na(coverpi(NULL, 0.5)))
})

# Coverage of confidence intervals
test_that("coverci returns expected output", {
  expect_equal(coverci(c(-0.5, 0.5)), 1)
  expect_equal(coverci(c(0.0, 0.5)), 0)
  expect_true(is.na(coverci(c(-0.5, NA))))
  expect_true(is.na(coverci(c(0))))
  expect_true(is.na(coverci(NULL)))
})

# Interval skewness
test_that("ski returns expected output", {
  expect_equal(ski(c(0, 2), 1), 0)
  expect_equal(ski(c(0, 3), 1), 0.3333333, tolerance = 1e-5)
  expect_true(is.na(ski(c(1, 1), 1)))  # zero width interval
  expect_true(is.na(ski(c(1, NA), 1)))
  expect_true(is.na(ski(c(1), 1)))
})

# Fishers skewness coefficient
test_that("fwskew computes skewness and handles edge cases", {
  es <- c(1.5, 2.3, 1.7)
  se <- c(0.11, 0.25, 0.16)
  expect_type(fwskew(es, se), "double")
  expect_true(is.na(fwskew(es, c(0, 1, 1))))  # division by zero
  expect_true(is.na(fwskew(es, c(NA, 0.1, 0.1))))
  expect_true(is.na(fwskew(es, c(0.1, 0.1))))
})

# Predictive probabilities
test_that("pg computes probability estimates correctly", {
  s <- c(NA, 0.1, 0.5, 0.9)
  expect_equal(pg(s, 0.5), 2 / 3)
  expect_equal(pg(s, 1), 0)
  expect_true(is.na(pg(NULL, 0.5)))
  expect_true(is.na(pg(c(NA, NA), 0.5)))
  expect_true(is.na(pg(s, NA)))
})

# CRPS
test_that("simcrps ignores NAs in both s and tn", {
  s <- c(NA, 1, 2)
  tn <- c(1, 2, NA)
  result <- simcrps(s, tn)
  expect_type(result, "double")
  expect_true(!is.na(result))
})

test_that("simcrps returns NA if all s or tn are NA", {
  expect_true(is.na(simcrps(c(NA, NA), c(1, 2))))
  expect_true(is.na(simcrps(c(1, 2), c(NA, NA))))
})

