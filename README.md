# sim-metaprediction
Simulation study investigating the predictive distributions for meta-analysis implemented in the `metaprediction` package. We provide a concise description based on the ADEMP framework.

## Aims
- Assess performance of predictive distributions in terms of calibration and sharpness.
- Assess performance of new estimator for the average effect.

## Data-Generating Process
We varied:
- the number of studies k ∈ {3, 5, 10, 20, 50}.
- the between-study heterogeneity determined by Higgins’ I2 ∈ {0.0, 0.3, 0.6, 0.9}.
- the number of large studies ∈ {0, 1, 2}.
- the underyling effect distribution ∈ {Normal, Skew-normal}.

Study sizes were set to 50 for normal studies and to 500 for large studies. We performed the simulation study in a full-factorial manner, investigating all combinations of factor levels. 

## Estimands and Other Targets
The true average effect was set to μ = −0.3.  Our aim to evaluate predictive performance represents a target, since it involves predicting the entire distribution of unobserved future effects rather than estimating a single parameter with a fixed value.

## Methods
Methods differ across performance measures:

- 95% prediction intervals: Equi-tailed 95% prediction intervals obtained from the three predictive distributions in the `metaprediction`package, the Higgins-Thompson-Spiegelhalter prediction interval and the parametric boostrap prediction  interval.
- Predictive distributions: The three predictive distributions in the `metaprediction`package and the Higgins-Thompson-Spiegelhalter predictive distribution.
- Point estimation of the average effect: Estimator for random-effects meta-analysis from the `metaprediction`package, the classical inverse-variance weights point estimator from fixed-effect and random-effects meta-analysis, and the point estimator from the Edgington combined p-value function without and with additive heterogeneity adjustment.
-  95% confidence intervals: Equi-tailed 95% confidence intervals for random-effects meta-analysis from the `metaprediction` package, classical random-effects meta-analysis, Hartung–Knapp method and Edgington’s method without and with additive heterogeneity adjustment.

## Performance Measures
- Coverage of 95% prediction intervals
- Skewness of 95% prediction intervals
- Width of 95% prediction intervals
- Continuous ranked probability scores
- Predictive probabilities
- Bias of the point estimator
- MSE of the point estimator
- Coverage of 95% confidence intervals
- Skewness of 95% confidence intervals
- Width of 95% confidence intervals
- Non-convergences
- Computing time

## Number of Iterations
For each of the 120 scenarios, 4000 simulation iterations were run to achieve a maximum Monte Carlo standard error of 0.005.

## Computational Details
The entire simulation study is programmed in the R programming language and conducted in R version 4.5.0 (2025-04-11) on a remote Debian GNU/Linux server (platform: x86_64-pc-linux-gnu). During the simulation, we performed the following steps iteratively:

- A set of factor levels is selected from a grid containing all possible combinations of factor levels.
-  The simulation iterations are distributed evenly across parallelized clusters. To ensure reproducibility of the simulation, we used random number generator (RNG) streams together with an intially set state of the RNG before initiation of the simulation. We use the doRNG package which uses the L’Ecuyer-CMRG RNG for the distribution. 
- After all parallel clusters have completed their iterations, the results are pooled and stored.
- The process is repeated for all sets of factor levels.
