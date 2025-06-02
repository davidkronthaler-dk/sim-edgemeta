# ------------------------------------------------------------------------------
# Utility functions for the simulation
# ------------------------------------------------------------------------------

## Catch errors
##------------------------------------------------------------------------------
trySim <- function(expr, nx = 1) {
  tryCatch(
    expr = expr,
    error = function(e)
      return(rep(NA_real_, nx))
  )
}

## Coverage of prediction interval
## -----------------------------------------------------------------------------
coverpi <- function(i, t) {
  c <- trySim({
    if (length(i) != 2 || anyNA(i) || is.null(i))
      return(NA_real_)
    mean(i[1] <= t & i[2] >= t, na.rm = TRUE)
  })
  return(c)
}

## Coverage of confidence interval
## -----------------------------------------------------------------------------
coverci <- function(i) {
  c <- trySim({
    if (length(i) != 2 || anyNA(i) || is.null(i))
      return(NA_real_)
    as.numeric(i[1] <= -0.3 & i[2] >= -0.3)
  })
  return(c)
}

## Skewness of interval
##------------------------------------------------------------------------------
ski <- function(i, p) {
  sk <- trySim({
    if (length(i) != 2 || anyNA(i) || i[2] == i[1])
      return(NA_real_)
    (i[2] + i[1] - 2 * p) / (i[2] - i[1])
  })
  return(sk)
}

## Fishers (weighted) skewness coefficient
##------------------------------------------------------------------------------

fwskew <- function(es, se) {
  fsk <- trySim({
    if (length(es) != length(se) ||
        anyNA(se) || anyNA(es) || any(se == 0))
      return(NA_real_)
    wi <- 1 / (se ^ 2)
    mu_bar <- sum(es * wi) / sum(wi)
    num <- (sum(wi * (es - mu_bar) ^ 3)) * sqrt(sum(wi))
    denum <- (sum(wi * (es - mu_bar) ^ 2)) ^ (3 / 2)
    if (denum == 0)
      return(NA_real_)
    num / denum
  })
  return(fsk)
}

## Predicted probability (larger equal to q)
##------------------------------------------------------------------------------

pg <- Vectorize(function(s, q) {
  trySim({
    if (is.null(s) || is.null(q) || !is.numeric(s) || all(is.na(s))) {
      return(NA_real_)
    }
    mean(s >= q, na.rm = TRUE)
  })
}, vectorize.args = "q")

## Continuous ranked probability score
##------------------------------------------------------------------------------

simcrps <- function(s, tn) {
  sc <- trySim({
    if (is.null(s) ||
        is.null(tn) || !is.numeric(s) || !is.numeric(tn))
      return(NA_real_)
    
    s_clean <- s[!is.na(s)]
    tn_clean <- tn[!is.na(tn)]
    
    if (length(s_clean) == 0 ||
        length(tn_clean) == 0)
      return(NA_real_)
    
    mean(sapply(as.list(tn_clean), function(x)
      scoringRules::crps_sample(x, s_clean)),
      na.rm = TRUE)
  })
  return(sc)
}
