# ------------------------------------------------------------------------------
# Functions computing performance measures for
# - predictive distributions (pd)
# - estimators for the average effect (mu)
# ------------------------------------------------------------------------------

performance_pd <- function(mm, tn, dt, k) {
  
  # 95% prediction intervals
  pis <- list(mm$pi.hts, mm$pi.nnf, mm$pi.fix, trySim({ mm$pd.simple$PI }, 2), 
              trySim({ mm$pd.full$PI },2))
  
  # Samples from pred. distribution of Higgins (2008)
  samphts <- trySim( {rt(1e5, k - 2) * sqrt(mm$p.meta1$tau2 + mm$p.meta1$seTE.random) + mm$p.meta1$TE.random}, 1e5)
  
  # Median of predictive distributions
  pdmed <- c(
    rep(trySim(mm$p.meta1$TE.random), 2),
    if (!is.null(mm$pd.fix$samples)) {
      trySim( {median(mm$pd.fix$samples[, "theta_new"], na.rm = TRUE) })
    } else trySim( {mm$pd.fix$estimate} ),
    trySim(median(mm$pd.simple$samples[,"theta_new"], na.rm = TRUE)),
    trySim(median(mm$pd.full$samples[,"theta_new"], na.rm = TRUE))
  )
  
  # Coverage of 95% prediction intervals
  picvr <- sapply(pis, function(x) {coverpi(x, tn)})
  
  # Width of 95% prediction intervals
  piwd <- sapply(pis, function(x) {
    if (is.numeric(x) && length(x) == 2 && all(!is.na(x))) x[2] - x[1] else NA_real_
  })
  
  # Skewness of 95% prediction intervals
  pisk <- mapply(function(x, p) {
    if (length(x) == 2 && !anyNA(x) && !is.na(p)) ski(x, p) else NA_real_
  }, x = pis, p = pdmed)
  
  # Predictive probabilities
  qhts <- trySim(pg(samphts, dt$q), 3)
  qfix <- if (is.list(mm$pd.fix) && !is.null(mm$pd.fix$samples)) {
    trySim( {pg(mm$pd.fix$samples[,"theta_new"], dt$q)} , 3)
  } else if (is.list(mm$pd.fix) && !is.null(mm$pd.fix$CI)) {
    rep(0, 3) # tau2 is estimated zero, no predictive distribution available
  } else {
    rep(NA_real_, 3)
  }
  qsimple <- trySim( {pg(mm$pd.simple$samples[,"theta_new"], dt$q)} , 3)
  qfull <- trySim({ pg(mm$pd.full$samples[,"theta_new"], dt$q) }, 3)
  
  # CRPS
  # 10'000 / 100'000 samples are used for computational reasons
  crpshts <- simcrps(sample(samphts, 1e4), tn)
  crpsfix <- if (is.list(mm$pd.fix) && !is.null(mm$pd.fix$samples)) {
    simcrps(sample(mm$pd.fix$samples[,"theta_new"], 1e4), tn)
  } else {
    simcrps(rep(mm$held2025a$estimate, 1e4), tn)
  }
  crpssimple <- simcrps(sample(mm$pd.simple$samples[,"theta_new"], 1e4), tn)
  crpsfull <- simcrps(sample(mm$pd.full$samples[, "theta_new"], 1e4), tn)
  
  # return
  return(list(
    pdmed = pdmed,
    picvr = picvr,
    piwd = piwd,
    pisk = pisk,
    qhts = qhts,
    qfix = qfix,
    qsimple = qsimple,
    qfull = qfull,
    crpshts = crpshts,
    crpsfix = crpsfix,
    crpssimple = crpssimple,
    crpsfull = crpsfull
  ))
}

performance_mu <- function(mm) {
  # Bias of point estimator for mu
  muhat <- c(trySim( {mm$p.meta1$TE.fixed} ),    
             trySim( {mm$p.meta1$TE.random} ),
             trySim( {mm$held2025u$estimate} ), 
             trySim( {mm$held2025a$estimate} ),
             trySim( {mean(mm$pd.full$samples[,"mu"], na.rm = T)} ))
  bias <- ifelse(is.na(muhat), NA_real_, muhat + 0.3)
  
  # Squared error of point estimator for mu
  sqe <- bias ^ 2
  
  # Coverage of 95% CI for mu
  cis <- list(trySim( {c(mm$p.meta2$lower.random, mm$p.meta2$upper.random)} , 2),
              trySim( {c(mm$p.meta1$lower.random, mm$p.meta1$upper.random)} , 2),
              trySim( {c(mm$held2025u$CI[1], mm$held2025u$CI[2])} , 2),
              trySim( {c(mm$held2025a$CI[1], mm$held2025a$CI[2])} , 2),
              trySim( {quantile(mm$pd.full$samples[, "mu"], p = c(0.025, 0.975), na.rm = T)} , 2))
  cicvr <- sapply(cis, function(x) {coverci(x)})
  
  # Width of 95% CI for mu
  ciwd <- sapply(cis, function(x) {
    if (is.numeric(x) && length(x) == 2 && all(!is.na(x))) x[2] - x[1] else NA_real_
  })
  
  # Skewness of 95% CI for mu
  cisk <- mapply(function(x, p) {
    if (length(x) == 2 && !anyNA(x) && !is.na(p)) ski(x, p) else NA_real_
  }, x = cis,
  p = c(trySim( {mm$p.meta1$TE.random} ),
        trySim( {mm$p.meta1$TE.random} ), 
        trySim( {mm$held2025u$estimate} ),
        trySim( {mm$held2025a$estimate} ),
        trySim( {mean(mm$pd.full$samples[,"mu"], na.rm = T)} )))
  
  # return
  return(list(
    bias = bias, 
    sqe = sqe,
    cicvr = cicvr,
    ciwd = ciwd,
    cisk = cisk
  ))
}