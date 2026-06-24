# sim-edgemeta

This repository contains code and results of a simulation study investigating Edgington's predictive distributions and the CD-Edgington estimator for random-effects meta-analysis. The investigated methods are implemented in the `edgemeta` package ([https://github.com/davidkronthaler-dk/edgemeta]). Further, the folder `manuscriptcode` contains the code for the manuscript:

Kronthaler, D., & Held, L. (2026). *Prediction intervals for random-effects meta-analysis based on confidence distributions and Edgington's method*. [arXiv. https://doi.org/10.48550/arXiv.2510.13216]
# Simulation Study
We provide a concise description based on the ADEMP framework.

## Aims
- Assess performance of Edgington's predictive distributions in terms of calibration and sharpness.
- Assess performance of CD-Edgington estimator for the average effect.

## Data-Generating Process
We varied:
- the number of studies k ∈ {3, 5, 10, 20, 50}.
- the degree of between-study heterogeneity ∈ {0%, 30%, 60%, 90%}.
- the number of large studies ∈ {0, 1, 2}.
- the underyling effect distribution ∈ {Normal, Skew-normal}.

Study sizes were set to 50 for normal studies and to 500 for large studies. We performed the simulation study in a full-factorial manner, investigating all combinations of factor levels. Effect estimates are generated on the standardized mean difference scale.

## Estimands and Other Targets
The true average effect was set to μ = −0.3.  Our aim to evaluate predictive performance represents a target, since it involves predicting the entire distribution of future effects rather than estimating a single parameter with a fixed value.

## Methods
Methods differ across performance measures:

- 95% prediction intervals: Equi-tailed 95% prediction intervals obtained from Edgington's predictive distributions in the `edgemeta`package, the Higgins-Thompson-Spiegelhalter prediction interval and the parametric boostrap prediction interval.
- Predictive distributions: Edgington's predictive distributions from the `edgemeta`package and the Higgins-Thompson-Spiegelhalter predictive distribution.
- 95% confidence intervals: CD-Edgington, classical random-effects meta-analysis, Hartung-Knapp and Edgington’s method without heterogeneity uncertainty adjustment.
- Point estimation of the average effect: CD-Edgington, classical inverse-variance weights estimator from random-effects meta-analysis, Edgington's method without heterogeneity uncertainty adjustment.

## Performance Measures
- Coverage of 95% confidence intervals
- Skewness of 95% confidence intervals
- Width of 95% confidence intervals
- Bias of the point estimator
- MSE of the point estimator
- Coverage of 95% prediction intervals
- Skewness of 95% prediction intervals
- Width of 95% prediction intervals
- Continuous ranked probability scores
- Non-convergences
- Computing time

## Number of Iterations
For each of the 120 scenarios, 4000 simulation iterations were run.

## Computational Details
The entire simulation study is programmed in the R programming language and conducted in R version 4.5.0 (2025-04-11) on a remote Debian GNU/Linux server (platform: x86_64-pc-linux-gnu). During the simulation, we performed the following steps iteratively:

- A set of factor levels is selected from a grid containing all possible combinations of factor levels.
- The simulation iterations are distributed evenly across parallelized clusters. To ensure reproducibility of the simulation, we used random number generator (RNG) streams together with an intially set state of the RNG before initiation of the simulation. We use the doRNG package which uses the L’Ecuyer-CMRG RNG for the distribution. 
- After all parallel clusters have completed their iterations, the results are pooled and stored.
- The process is repeated for all sets of factor levels.
