#-------------------------------------------------------------------------------
# Data-generating function
# Function for one simulation iteration
# Function for multiple iterations
# Function for parallelization
# Function for iterating over scenarios
#-------------------------------------------------------------------------------

## Data-generating process 
##------------------------------------------------------------------------------
dgp <- function(k, I2, k_large, dist) {
  
  # Study sizes
  ni <- c(rep(500, k_large), rep(50, k - k_large))
  
  # Average effect
  mu <- -0.3
  
  # Simulate standard errors
  se <- sqrt(stats::rchisq(k, df = 2 * (ni - 1)) * ((ni - 1) * ni)^(-1))
  
  # Compute tau2
  tau2 <- 1 / k * sum(2 / ni) * (I2 / (1 - I2))
  
  # Simulate true effects and compute quantiles from effect distribution
  if (dist  == "LSN") {
    alpha <- - 4
    d <- alpha / sqrt(1 + alpha ^ 2)
    omega <- sqrt(tau2 / (1 - 2 * (d ^ 2) / pi))
    xi <- mu - omega * d * sqrt(2 / pi)
    es <- sn::rsn(k + 10000, xi = xi, omega = omega, alpha = alpha)
    q <- sn::qsn(c(0.5, 0.9, 0.95), xi = xi, omega = omega + .Machine$double.eps, alpha = alpha)
    
  } else if (dist == "N") {
    es <- stats::rnorm(k + 10000, mean = mu, sd = sqrt(tau2))
    q <- stats::qnorm(c(0.5, 0.9, 0.95), mean = mu, sd = sqrt(tau2))
  }
  
  # Simulate effect estimates
  hes <- rnorm(k, es, sqrt(2/ni))
  
  # Return
  return(list(hes = hes, se = se, es = es, q = q))
  
}

## One Simulation Iteration
##------------------------------------------------------------------------------

onesim <- function(k, I2, k_large, dist) {
  
  # Generate data 
  dt <- dgp(k = k, I2 = I2, k_large = k_large, dist = dist)
  
  # Extract future effects
  tn <- trySim( {dt$es[(k + 1):length(dt$es)]}, 10000)
  
  # Fishers weighted skewness for effect estimates
  skhes <- fwskew(dt$hes, dt$se)
  
  # Fishers skewness for effects
  skes <- if (I2 == 0) { 0 } else { fwskew(dt$es[1:k], rep(1, k)) }
  
  # Vector to store computation time
  t <- rep(NA, 5)
  
  # ----- Fit methods ----- #
  # 'meta' package
  t[1] <- system.time(p.meta1 <- trySim({ metagen(TE = dt$hes, seTE = dt$se, 
                                                  random = T, method.tau2 = "REML", 
                                                  method.random.ci = "HK",
                                                  method.predict = "HTS") }
                                        ))["elapsed"]
  
  t[2] <- system.time(p.meta2 <- trySim({ metagen(TE = dt$hes, seTE = dt$se, 
                                                  random = T, method.tau2 = "REML", 
                                                  method.random.ci = "HK",
                                                  method.predict = "NNF", 
                                                  B = 100000) }
                                        ))["elapsed"]
  
  # Store HTS and NNF prediction interval
  pi.hts <- trySim({if (!is.null(p.meta1$lower.predict) && !is.null(p.meta1$upper.predict)) 
    c(p.meta1$lower.predict, p.meta1$upper.predict) else rep(NA_real_, 2)}, 2)
  
  pi.nnf <- trySim({if (!is.null(p.meta2$lower.predict) && !is.null(p.meta2$upper.predict)) 
    c(p.meta2$lower.predict, p.meta2$upper.predict) else rep(NA_real_, 2)}, 2)
  
  # 'confMeta' package (Held et al., 2025)
  held2025 <- trySim({ confMeta(estimates = dt$hes, SEs = dt$se, 
                                fun = p_edgington,
                                fun_name = "Edgington  (one-sided input)",
                                input_p = "one.sided") })
  
  # 'metaprediction' package
  t[3] <- system.time(pd.fix <- trySim({ PredDist(es = dt$hes, se = dt$se, 
                                                  method = "FixedTau2") }
                                       ))["elapsed"]
  
  pi.fix <- trySim({ if (!is.null(pd.fix$PI)) {pd.fix$PI} else {pd.fix$CI} }, 2)
  
  t[4] <- system.time(pd.simple <- trySim({ PredDist(es = dt$hes, se = dt$se, 
                                                     method = "SimplifiedCD") }
                                          ))["elapsed"]
  
  t[5] <- system.time(pd.full <- trySim({ PredDist(es = dt$hes, se = dt$se, 
                                                   method = "FullCD") }
                                        ))["elapsed"]
  
  # 95% prediction intervals
  pis <- list(pi.hts, pi.nnf, pi.fix, trySim({ pd.simple$PI }, 2), trySim({ pd.full$PI },2))
  
  # Samples from pred. distribution of Higgins (2008)
  samphts <- trySim( {rt(1e5, k - 2) * sqrt(p.meta1$tau2 + p.meta1$seTE.random) + p.meta1$TE.random})
  
  # ----- Analysis ----- #
  # Median of predictive distributions
  pdmed <- c(
    rep(trySim(p.meta1$TE.random), 2),
    trySim(held2025$p_max[,"x"]),
    trySim(median(pd.simple$samples[,"theta_new"], na.rm = TRUE)),
    trySim(median(pd.full$samples[,"theta_new"], na.rm = TRUE))
  )
  
  # Coverage of 95% prediction intervals
  picvr <- sapply(pis, function(x) {coverpi(x, tn)})
  
  # Width of 95% prediction intervals
  piwd <- sapply(pis, function(x) {
    if (is.numeric(x) && length(x) == 2 && all(!is.na(x))) x[2] - x[1] else NA_real_
  })
  
  # Skewness of 95% prediction intervals
  pisk <- pisk <- mapply(function(x, p) {
    if (length(x) == 2 && !anyNA(x) && !is.na(p)) ski(x, p) else NA_real_
  }, x = pis, p = pdmed)
  
  # Predictive probabilities
  qhts <- trySim(pg(samphts, dt$q), 3)
  qfix <- if (is.list(pd.fix) && !is.null(pd.fix$samples)) {
    trySim( {pg(pd.fix$samples[,"theta_new"], dt$q)} , 3)
  } else if (is.list(pd.fix) && !is.null(pd.fix$CI)) {
    rep(0, 3) # tau2 is estimated zero, no predictive distribution available
  } else {
    rep(NA_real_, 3)
  }
  qsimple <- trySim( {pg(pd.simple$samples[,"theta_new"], dt$q)} , 3)
  qfull <- trySim({ pg(pd.full$samples[,"theta_new"], dt$q) }, 3)
  
  # CRPS
  # 10'000 / 100'000 samples are used for computational reasons
  crpshts <- simcrps(sample(samphts, 1e4), tn)
  crpsfix <- if (is.list(pd.fix) && !is.null(pd.fix$samples)) {
    simcrps(sample(pd.fix$samples[,"theta_new"], 1e4), tn)
  } else {
    simcrps(held2025$p_max[,"x"], tn)
  }
  crpssimple <- simcrps(sample(pd.simple$samples[,"theta_new"], 1e4), tn)
  crpsfull <- simcrps(sample(pd.full$samples[, "theta_new"], 1e4), tn)

  # Bias of point estimator for mu
  muhat <- c(trySim( {p.meta1$TE.random} ),
            trySim( {held2025$p_max[,"x"]} ), 
            trySim( {mean(pd.full$samples[,"mu"], na.rm = T)} )) 
  bias <- ifelse(is.na(muhat), NA_real_, muhat + 0.3)
  
  # Squared error of point estimator for mu
  sqe <- bias ^ 2
  
  # Coverage of 95% CI for mu
  cis <- list(trySim( {c(p.meta1$lower.random, p.meta1$upper.random)} , 2),
              trySim( {held2025$joint_cis} , 2),
              trySim( {quantile(pd.full$samples[, "mu"], p = c(0.025, 0.975), na.rm = T)} , 2))
  cicvr <- sapply(cis, function(x) {coverci(x)})
  
  # Width of 95% CI for mu
  ciwd <- sapply(cis, function(x) {
    if (is.numeric(x) && length(x) == 2 && all(!is.na(x))) x[2] - x[1] else NA_real_
  })
  
  # Skewness of 95% CI for mu
  cisk <- mapply(function(x, p) {
    if (length(x) == 2 && !anyNA(x) && !is.na(p)) ski(x, p) else NA_real_
  }, x = cis,
  p = c(trySim( {p.meta1$TE.random} ), 
        trySim( {held2025$p_max[, "x"]} ),
        trySim( {mean(pd.full$samples[,"mu"], na.rm = T)} )))
  
  # ----- Return -----
  r <- c(picvr, piwd, pisk, t, # Coverage, width, skewness, comp. time of 5 prediction interval
         qhts, qfix, qsimple, qfull, # Predicted probabilities 4 predictive dist.
         crpshts, crpsfix, crpssimple, crpsfull, # CRPS 4 predictive dist.
         bias, sqe, # Bias and squared error for 3 point estimators,
         cicvr, ciwd, cisk, # Coverage, width and skewness of CIs for 3 methods
         skhes, skes, trySim( {p.meta1$tau2} )) # Skewness effect estimates, effects, estimated tau2
  
  # Return
  return(r)
}

## Function for simulation study
##------------------------------------------------------------------------------

simstudy <-  function(k, I2, k_large, dist, niter) {

  # Perform iterations
  res_list <- replicate(niter, trySim({ onesim(
    k = k,
    I2 = I2,
    dist = dist,
    k_large = k_large
  )}, 54), simplify = FALSE)

  # Transform to matrix
  res <- tryCatch(
    do.call(rbind, res_list),
    error = function(e) matrix(NA_real_, nrow = niter, ncol = 54)
  )

  # Factor levels
  fl <- data.frame(k = k, I2 = I2, dist = dist, k_large = k_large)
  
  # Return object
  r <- cbind(fl, res)
  colnames(r)[5:58] <- 5:58

  # Return
  return(r)
}

## Function for parallelization
##------------------------------------------------------------------------------
parsimstudy <- function(k, I2, k_large, dist, niter) {
  
  # Number of parallelized CPUs
  n_cores <- 25 
  
  # Number of iterations per core
  niter_per_cpu <- ceiling(niter/n_cores)
  
  # Construct cluster
  cl = parallel::makeCluster(n_cores)
  
  # Shutdown the cluster on exit
  on.exit(parallel::stopCluster(cl))
  
  # Export functions to cluster
  clusterExport(cl, varlist = c("simstudy", "onesim", "dgp", "coverpi", "coverci", 
                                "ski", "fwskew", "pg", "simcrps", "trySim"))
  
  # Register for parallel processing
  registerDoParallel(cl)
  
  # Parallelized computation
  res <- foreach(core_nr = 1:n_cores,
                 .packages = c("meta", "sn", "metaprediction", "confMeta", "scoringRules"),
                 .combine = rbind) %dorng% {
                   simstudy(k = k, 
                            I2 = I2, 
                            k_large = k_large, 
                            dist = dist,
                            niter = niter_per_cpu)
                 }
  
  # Return
  return(res)
}

## Function to iterate over parameter grid
##------------------------------------------------------------------------------

itersim <- function(k, I2, k_large, dist, niter) {
  
  # Factor levels
  paras <- expand.grid(
    k = k,
    I2 = I2,
    dist = dist,
    k_large = k_large
  )
  
  # Iterate over parameter combinations
  for (param in 1:nrow(paras)) {
    
    # Console message
    message(
      "\n[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
      "Running Parameter Iteration ", param, " of ", nrow(paras), "\n",
      "Parameters: ",
      paste(names(paras[param, ]), paras[param, ], sep = " = ", collapse = ", ")
    )
    
    # Catch errors
    try({
      res <- parsimstudy(k = paras$k[param],
                         I2 = paras$I2[param],
                         dist = paras$dist[param],
                         k_large = paras$k_large[param],
                         niter = niter)
      
      colnames(res)[5:58] <- c("pi.cvr.hts", "pi.cvr.nnf", "pi.cvr.fix", "pi.cvr.simple", "pi.cvr.full",
                               "pi.w.hts", "pi.w.nnf", "pi.w.fix", "pi.w.simple", "pi.w.full",
                               "pi.sk.hts", "pi.sk.nnf", "pi.sk.fix", "pi.sk.simple", "pi.sk.full",
                               "t.hts", "t.nnf", "t.fix", "t.simple", "t.full",
                               "q05.hts", "q01.hts", "q005.hts",
                               "q05.fix", "q01.fix", "q005.fix",
                               "q05.simple", "q01.simple", "q005.simple",
                               "q05.full", "q01.full", "q005.full",
                               "crps.hts", "crps.fix", "crps.simple", "crps.full",
                               "bias.ivw", "bias.held", "bias.samp",
                               "sqe.ivw", "sqe.held", "sqe.samp",
                               "ci.cvr.hk", "ci.cvr.held", "ci.cvr.samp",
                               "ci.w.hk", "ci.w.held", "ci.w.samp",
                               "ci.sk.hk", "ci.sk.held", "ci.sk.samp", 
                               "sk.hes", "sk.es", "es.tau2")
                         
      saveRDS(res, file = paste0("../results/sr",
                                 paste(gsub("\\.", "", paras[param, ]), collapse = "_"),
                                 ".RDS"))
    })
  }
}



























