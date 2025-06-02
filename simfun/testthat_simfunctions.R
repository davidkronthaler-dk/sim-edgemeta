#-------------------------------------------------------------------------------
# Unit tests for simulation functions
#-------------------------------------------------------------------------------

# Data-generating mechanism
test_that("data-generating mechanisms works for all scenarios", {
  fltest <- expand.grid(
    k = c(3, 5, 10, 20, 50),
    I2 = c(0, 0.3, 0.6, 0.9),
    k_large = c(0, 1, 2),
    dist = c("N", "LSN")
  )
  
  apply(fltest, 1, function(row) {
    expect_no_error(dgp(
      k = as.numeric(row["k"]),
      I2 = as.numeric(row["I2"]),
      k_large = as.numeric(row["k_large"]),
      dist = row["dist"]
    ))
  })
  
})

# Simstudy (iterations for one scenario)
test_that("simstudy doesn't break with invalid niter", {
  result <- suppressWarnings(simstudy(0, 0, 0, "N", 1))
  
  # Check that a data frame is returned
  expect_s3_class(result, "data.frame")
  
  # Should have 1 row (since replicate fails and fallback matrix is used)
  expect_equal(nrow(result), 1)
  
  # Should have 4 metadata columns + 54 NA columns
  expect_equal(ncol(result), 58)
  
})
