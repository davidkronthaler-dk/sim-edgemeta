# ------------------------------------------------------------------------------
# Main document for the simulation study on random-effect meta-analysis evaluating 
# - Predictive distributions
# - 95% Prediction intervals 
# - Point estimators for the average effect
# - Confidence intervals for the average effect
# We compare methods proposed by us to standard methods implemented in the 
# 'meta' package and methods proposed by Held et al. (2025) 
#-------------------------------------------------------------------------------

## Load libraries
##------------------------------------------------------------------------------
library(doRNG)             # RNG streams
library(doParallel)        # Parallelization
library(metaprediction)    # Our package
library(sn)                # Skew-normal distribution
library(meta)              # Standard methods
library(confMeta)          # Methods by Held et al. (2025)
library(scoringRules)      # Continuous ranked probability scores
library(testthat)          # Tests for simulation functions

## Source functions and perform unit tests
##------------------------------------------------------------------------------
source("simutils.R")
source("testthat_utils.R")
source("simfunctions.R")
source("testthat_simfunctions.R")

## Run the simulation
##------------------------------------------------------------------------------
set.seed(021098)
itersim(k = c(3, 5, 10, 20, 50), I2 = c(0.0, 0.3, 0.6, 0.9),
        k_large = c(0, 1, 2), dist = c("N", "LSN"), niter = 4000)

## Store the computational details
## -----------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "../compinfo/sessionInfo.txt")

