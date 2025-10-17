# ------------------------------------------------------------------------------
# Functions to fit methods
# ------------------------------------------------------------------------------

metamethods <- function(dt) {
  
  # Vector to store computation time
  t <- rep(NA, 5)
  
  # 'meta' package
  t[1] <- system.time(p.meta1 <- trySim({ metaclassic(dt$hes, dt$se, pi = "HTS", cit = "HK") }
  ))["elapsed"]
  
  t[2] <- system.time(p.meta2 <- trySim({ metaclassic(dt$hes, dt$se, pi = "NNF", cit = "classic") }
  ))["elapsed"]
  
  # Store HTS and NNF prediction interval
  pi.hts <- trySim({if (!is.null(p.meta1$lower.predict) && !is.null(p.meta1$upper.predict)) 
    c(p.meta1$lower.predict, p.meta1$upper.predict) else rep(NA_real_, 2)}, 2)
  
  pi.nnf <- trySim({if (!is.null(p.meta2$lower.predict) && !is.null(p.meta2$upper.predict)) 
    c(p.meta2$lower.predict, p.meta2$upper.predict) else rep(NA_real_, 2)}, 2)
  
  # Methods from Held et al., 2025
  held2025u <- trySim({ remaeffect(dt$hes, dt$se, "NHEU") }) # 'held2025u' deprecated, not used
  held2025a <- trySim({ remaeffect(dt$hes, dt$se, "NHEU") })
  
  # 'edgemeta' package
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
  
  return(list(t = t, p.meta1 = p.meta1, p.meta2 = p.meta2, pi.hts = pi.hts,
              pi.nnf = pi.nnf, held2025u = held2025u, held2025a = held2025a,
              pd.fix = pd.fix, pi.fix = pi.fix, pd.simple = pd.simple, 
              pd.full = pd.full))
}


metaclassic <- function(es, se, pi, cit = "HK") {
  # Attempt 1: default settings
  result <- tryCatch(
    meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = "REML",
                  method.random.ci = cit, method.predict = pi, B = 100000
    ),
    error = function(e) NULL
  )
  
  # Attempt 2: increase max iterations
  if (is.null(result)) {
    result <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = "REML",
                    method.random.ci = cit, method.predict = pi, B = 100000,
                    control = list(maxiter = 10000)),
      error = function(e) NULL
    )
  }
  
  # Attempt 3: decrease step size
  if (is.null(result)) {
    result <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = "REML",
                    method.random.ci = cit, method.predict = pi, B = 100000,
                    control = list(stepadj = 0.5)),
      error = function(e) NULL
    )
  }
  
  # Attempt 4: increase both maxiter and decrease step size
  if (is.null(result)) {
    result <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = "REML",
                    method.random.ci = cit, method.predict = pi, B = 100000,
                    control = list(maxiter = 10000, stepadj = 0.5)),
      error = function(e) NULL
    )
  }
  
  # Attempt 5: increase both maxiter and decrease step size further
  if (is.null(result)) {
    result <- tryCatch(
      meta::metagen(TE = es, seTE = se, random = TRUE, method.tau = "REML",
                    method.random.ci = cit, method.predict = pi, B = 100000,
                    control = list(maxiter = 100000, stepadj = 0.25)),
      error = function(e) NULL
    )
  }
  
  return(result)
}
