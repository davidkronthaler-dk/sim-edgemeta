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

