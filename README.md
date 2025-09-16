# sim-edgemeta
Simulation study investigating the performance of the CD-Edgington estimator for the average effect and predictive distributions method for random-effects meta-analysis implemented in the `edgemeta` package. We provide a concise description based on the ADEMP framework.

## Aims
- Assess performance of CD-Edgington estimator for the average effect.
- Assess performance of predictive distributions in terms of calibration and sharpness.

## Data-Generating Process
We varied:
- the number of studies k ∈ {3, 5, 10, 20, 50}.
- the between-study heterogeneity determined by Higgins’ I2 ∈ {0.0, 0.3, 0.6, 0.9}.
- the number of large studies ∈ {0, 1, 2}.
- the underyling effect distribution ∈ {Normal, Skew-normal}.

Study sizes were set to 50 for normal studies and to 500 for large studies. We performed the simulation study in a full-factorial manner, investigating all combinations of factor levels. Effect estimates are generated on the standardized mean difference scale.

## Estimands and Other Targets
The true average effect was set to μ = −0.3.  Our aim to evaluate predictive performance represents a target, since it involves predicting the entire distribution of unobserved future effects rather than estimating a single parameter with a fixed value.

## Methods
Methods differ across performance measures:

-  95% confidence intervals: CD-Edgington, classical random-effects meta-analysis, Hartung–Knapp and Edgington’s method with additive heterogeneity adjustment.
- Point estimation of the average effect: CD-Edgington, classical inverse-variance weights estimator from random-effects meta-analysis, Edgington's method with additive heterogeneity adjustment.
- 95% prediction intervals: Equi-tailed 95% prediction intervals obtained from the three predictive distributions in the `edgemeta`package, the Higgins-Thompson-Spiegelhalter prediction interval and the parametric boostrap prediction  interval.
- Predictive distributions: The three predictive distributions in the `edgemeta`package and the Higgins-Thompson-Spiegelhalter predictive distribution.

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
-  The simulation iterations are distributed evenly across parallelized clusters. To ensure reproducibility of the simulation, we used random number generator (RNG) streams together with an intially set state of the RNG before initiation of the simulation. We use the doRNG package which uses the L’Ecuyer-CMRG RNG for the distribution. 
- After all parallel clusters have completed their iterations, the results are pooled and stored.
- The process is repeated for all sets of factor levels.
