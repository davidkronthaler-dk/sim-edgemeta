## Clear environment
##------------------------------------------------------------------------------
rm(list = ls())

## Library Packages
## -----------------------------------------------------------------------------
# the following packages are available from CRAN 
# (install with install.packages("PACKAGE"))
library(dplyr) 
library(tidyr)
library(ggplot2)
library(ggthemes)
library(latex2exp)
library(sn)
library(xtable)
library(patchwork)
library(meta)
library(confMeta)
library(metafor)
library(ggtext)
library(numDeriv)
library(parallel)
library(pimeta)
library(scales)

# the following packages are available from Github
# (install with
#  remotes::install_github(repo = "EBPI-Biostatistics/biostatUZH", build_vignettes = TRUE)
# )
library(biostatUZH)

# (install with 
#  remotes::install_github("davidkronthaler-dk/edgemeta)
# )
library(edgemeta)

## Additional settings
## -----------------------------------------------------------------------------
options(width = 85, digits = 4, show.signif.stars = FALSE)

## Source R Code
##------------------------------------------------------------------------------
source("functions/utility_functions.R")

## Plotting
##------------------------------------------------------------------------------
# ggplot2 theme
theme_dk <- function() {
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linetype = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(colour = "black", face = "plain",
                             family = "Times",
                             size = 10),
    axis.title = element_text(colour = "black",
                              size = 10, face = "plain"),
    axis.ticks = element_line(colour = "black"),
    title = element_text(size = 10, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = "gray90", color = "black"), 
  ) 
}

# Axis labels
lab_mu <- expression(mu) 
lab_fmu <- expression(c(mu)) 

## Load data Covid example and Serenoa example
##------------------------------------------------------------------------------
load("data/Covid_example.Rdata")
load("data/Serenoa_example.Rdata")


## Web Appendix A: P-Value Functions and Confidence Distributions: Example 
## and Applications
## -----------------------------------------------------------------------------
# Estimate and grid for parameter mu
hmu <- data_covid_example$logOR[1]
sehmu <- data_covid_example$logSE[1]
h0 <- seq(-4.5, 6, l = 250)

# 95% confidence interval
upex1 <- hmu + qnorm(0.975) * sehmu
lpex1 <- hmu - qnorm(0.975) * sehmu

# One-sided p-value function for alternative "greater than"
seg_onesided <- data.frame(
  x = c(hmu, -Inf, upex1, -Inf, lpex1, -Inf),
  xend = c(hmu, hmu, upex1, upex1, lpex1, lpex1),
  y = c(-Inf, 0.5, 0, 0.975, 0, 0.025),
  yend = c(0.5, 0.5, 0.975, 0.975, 0.025, 0.025),
  line_type = factor(c("Median estimate", "Median estimate", 
                       "95% Confidence interval", "95% Confidence interval", 
                       "95% Confidence interval", "95% Confidence interval"),
                     levels = c("Median estimate", "95% Confidence interval"))
)

ponesided <- data.frame(mu = h0, p = p_fun_wald(h0, hmu, sehmu)) |> 
  ggplot(aes(x = mu, y = p)) +
  geom_vline(xintercept = c(-2, 0, 2), color = "gray90") +
  geom_hline(yintercept = c(0, 0.2, 0.4, 0.6, 0.8, 1), color = "gray90") +
  geom_line() +
  labs(y = expression(p["1s,+"](mu) == C(mu)), x = lab_mu, title = "A") +
  theme_dk() +
  theme(legend.text = element_text(size = 10)) +
  scale_y_continuous(
    expand = c(0.01, 0),
    breaks = c(0.025, 0.2, 0.4, 0.5, 0.6, 0.8, 0.975),
    labels = c("<span style='color:#3474b4;'>0.025</span>", "0.2", "0.4", 
               "<span style='color:#3474b4;'>0.500</span>", "0.6", "0.8",
               "<span style='color:#3474b4;'>0.975</span>")) +
  scale_x_continuous(
    expand = c(0.0, 0),
    breaks = c(-4, -2, lpex1,0, hmu, 2, upex1, 4),
    labels = c(
      "-4", "-2",
      paste0("<br><span style='color:#3474b4;'>", 
             round(lpex1, digits = 2), "</span>"),
      "0",
      paste0("<br><span style='color:#3474b4;'>", 
             round(hmu, digits = 2), "</span>"),
      "2",
      paste0("<br><span style='color:#3474b4;'>", 
             round(upex1, digits = 2), "</span>"),
      "4"
    )) +
  geom_segment(data = seg_onesided, 
               aes(x = x, xend = xend, y = y, yend = yend, linetype = line_type),
               col = "#3474b4") +
  scale_linetype_manual(values = c("Median estimate" = 1, 
                                   "95% Confidence interval" = 2)) +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown(),
        panel.grid.major = element_blank())


# Two-sided p-value function
seg_ptwosided <- data.frame(
  x = c(hmu, -Inf, lpex1, upex1, upex1, hmu),
  xend = c(hmu, upex1, lpex1, upex1, upex1, -Inf),
  y = c(-Inf, 0.05, -Inf, -Inf, 0.05, 1),
  yend = c(1, 0.05, 0.05, 0.05, 0.05, 1),
  line_type = factor(
    c("Median estimate",
      "95% Confidence interval", "95% Confidence interval",
      "95% Confidence interval", "95% Confidence interval",
      "Median estimate"),
    levels = c("Median estimate", "95% Confidence interval")
  )
)
  
ptwosided <- data.frame(mu = h0, p = p_fun_wald(h0, hmu, sehmu, F)) |>
  ggplot(aes(x = mu, y = p)) +
  geom_vline(xintercept = c(-2, 0, 2), color = "gray90") +
  geom_hline(yintercept = c(0, 0.2, 0.4, 0.6, 0.8, 1), color = "gray90") +
  geom_line() +
  labs(y = expression(p["2s"](mu)), x = lab_mu, title = "B") +
  theme_dk() +
  theme(legend.text = element_text(size = 10)) +
  scale_y_continuous(
  expand = c(0.01, 0),
  breaks = c(0.05, 0.2, 0.4, 0.6, 0.8, 1),
  labels = c("<span style='color:#3474b4;'>0.05</span>", "0.2", "0.4", 
             "0.6", "0.8", "<span style='color:#3474b4;'>1</span>")) +
  scale_x_continuous(
    expand = c(0.0, 0),
    breaks = c(-4, -2, lpex1,0, hmu, 2, upex1, 4),
    labels = c(
      "-4", "-2",
      paste0("<br><span style='color:#3474b4;'>", 
             round(lpex1, digits = 2), "</span>"),
      "0",
      paste0("<br><span style='color:#3474b4;'>", 
             round(hmu, digits = 2), "</span>"),
      "2",
      paste0("<br><span style='color:#3474b4;'>", 
             round(upex1, digits = 2), "</span>"),
      "4"
    )) +
   geom_segment(
    data = seg_ptwosided,
    aes(x = x, xend = xend, y = y, yend = yend, linetype = line_type),
    color = "#3474b4"
  ) +
  scale_linetype_manual(values = c("Median estimate" = 1, 
                                   "95% Confidence interval" = 2)) +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown(),
        panel.grid.major = element_blank())


# Confidence density
pcd <- data.frame(
  h0 = h0, 
  cd = numDeriv::grad(function(x) p_fun_wald(x, hmu, sehmu), h0)
) |>
  ggplot(aes(x = h0, y = cd)) + 
  geom_line() +
  labs(x = lab_mu, y = lab_fmu, title = "C") +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dk()

# Align plots
(ponesided + ptwosided + pcd) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(legend.position = "bottom"))


## Web Appendix B: The CD-Edgington estimator
## B.3 Illustrative Examples of Heterogeneity Confidence Densities
## Covid-19 example
## -----------------------------------------------------------------------------
covid_ma_meta <- metagen(data_covid_example$logOR, data_covid_example$logSE, 
                         method.tau = "PM")

# Confidence distribution function 
sxi_ctau2_covid <- seq(0, 2, l = 1000)
Ctau2_covid <- sapply(sxi_ctau2_covid, function(xi) {
  result <- integrate( function(x) {
    cdtau2(x, data_covid_example$logOR, data_covid_example$logSE)},
    lower = 0, upper = xi,
    rel.tol = 1e-8, abs.tol = 1e-12)
  result$value
})
Ctau2_covid <- Ctau2_covid / max(Ctau2_covid)

# Median and 95% confidence interval
QuantTau2_covid <- approxfun(Ctau2_covid, sxi_ctau2_covid, rule = 2)
median_q_covid <- QuantTau2_covid(0.5)
q025_covid <- QuantTau2_covid(0.025)
q975_covid <- QuantTau2_covid(0.975)

# Plot confidence density
p_ctau2_covid <- data.frame(
  xi = sxi_ctau2_covid,
  cd = cdtau2(sxi_ctau2_covid, data_covid_example$logOR, data_covid_example$logSE)) |>
  ggplot(aes(x = xi, y = cd)) +
  geom_line() +
  geom_segment(aes(x = median_q_covid, xend = median_q_covid, y = -Inf, yend = Inf,
               linetype = "Median", color = "Median"), linewidth = 0.3,
               data = data.frame()) +
  geom_segment(aes(x = q025_covid, xend = q025_covid, y = -Inf, yend = Inf,
                   linetype = "95% CI", color = "95% CI"), linewidth = 0.3,
               data = data.frame()) +
  geom_segment(aes(x = q975_covid, xend = q975_covid, y = -Inf, yend = Inf,
                   linetype = "95% CI", color = "95% CI"), linewidth = 0.3,
               data = data.frame()) +
  scale_y_continuous(expand = c(0.005,0)) +
  scale_x_continuous(
    expand = c(0.005, 0),
    breaks = c(q025_covid, median_q_covid, 0.5, 1, 1.5, q975_covid, 2),
    labels = c(
      paste0("<span style='color:#3474b4;'>", format(q025_covid, digits = 2),
             "</span>"),
      paste0("<span style='color:#3474b4;'><br>",
             format(median_q_covid, digits = 2),
             "</span>"),
      "0.5", "1", "1.5",
      paste0("<span style='color:#3474b4;'><br>",
             format(q975_covid, digits = 3),
             "</span>"),
      "2")) +
  scale_linetype_manual(
    name = NULL,
    values = c("Median" = 1, "95% CI" = 2),
    breaks = c("Median", "95% CI")
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Median" = "#3474b4", "95% CI" = "#3474b4"),
    breaks = c("Median", "95% CI")
  ) +
  labs(
    y = expression(c~(tau^2)), x = expression(tau^2),
    title = "A"
  ) + 
  theme_dk() +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown()) 

## Serenoa example
## -----------------------------------------------------------------------------
sere_ma_meta  <- metagen(serenoa$effect.size, serenoa$se, method.tau = "PM")

# Confidence distribution function
sxi_ctau2_sere <- seq(0, 5, l = 1000)
Ctau2_sere <- sapply(sxi_ctau2_sere, function(xi) {
  result <- integrate( function(x) {
    cdtau2(x, serenoa$effect.size, serenoa$se)},
    lower = 0, upper = xi,
    rel.tol = 1e-8, abs.tol = 1e-12)
  result$value
})
Ctau2_sere <- Ctau2_sere / max(Ctau2_sere)

# Median and 95% confidence interval
QuantTau2_sere <- approxfun(Ctau2_sere, sxi_ctau2_sere, rule = 2)
median_q <- QuantTau2_sere(0.5)
q025 <- QuantTau2_sere(0.025)
q975 <- QuantTau2_sere(0.975)

# Compute and plot confidence density
p_ctau2_sere <- data.frame(
  xi = sxi_ctau2_sere,
  cd = cdtau2(sxi_ctau2_sere, serenoa$effect.size, serenoa$se)) |>
  ggplot(aes(x = xi, y = cd)) +
  geom_segment(aes(x = median_q, xend = median_q, y = -Inf, yend = Inf,
                   linetype = "Median", color = "Median"),
               linewidth = 0.3, data = data.frame()) +
  geom_segment(aes(x = q025, xend = q025, y = -Inf, yend = Inf,
                   linetype = "95% CI", color = "95% CI"),
               linewidth = 0.3, data = data.frame()) +
  geom_segment(aes(x = q975, xend = q975, y = -Inf, yend = Inf,
                   linetype = "95% CI", color = "95% CI"),
               linewidth = 0.3, data = data.frame()) +
  geom_line() +
  scale_y_continuous(limits = c(0,0.9), expand = c(0,0)) +
  scale_x_continuous(
    expand = c(0.005, 0),
    breaks = c(0, q025, median_q, 1, 2, 3, q975, 4, 5),
    labels = c("0",
               paste0(
                 "<span style='color:#3474b4;'><br>",
                 format(q025, digits = 2), "</span>"),
               paste0("<span style='color:#3474b4;'><br>",
                      format(median_q, digits = 2), "</span>"),
               "1", "2", "3",
               paste0("<span style='color:#3474b4;'><br>",
                      format(q975, digits = 3), "</span>"),
               "4", "5")) +
  labs(
    y = expression(c~(tau^2)), x = expression(tau^2),
    title = "B"
  ) + 
  scale_linetype_manual(
    name = NULL,
    values = c("Median" = 1, "95% CI" = 2),
    breaks = c("Median", "95% CI")
  ) +
  scale_color_manual(
    name = NULL,
    values = c("Median" = "#3474b4", "95% CI" = "#3474b4"),
    breaks = c("Median", "95% CI")
  ) +
  theme_dk() +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown())

# Align plots
(
  p_ctau2_covid + p_ctau2_sere +
    plot_layout(axis_titles = "collect_x", guides = "collect")
) &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  )


## B.4 Comparison of Monte Carlo and GAQ Integration
## -----------------------------------------------------------------------------
if (!file.exists("data/MCvsGAQ.RDS")) {
  
  library(doParallel)
  library(foreach)
  
  # One simulation iteration
  one_simulation <- function(k, I2) {
    
    ni <- 50
    mu <- -0.3
    se <- sqrt(rchisq(k, df = 2 * (ni - 1)) * ((ni - 1) * ni)^(-1))
    tau2 <- 1 / k * sum(2 / ni) * (I2 / (1 - I2))
    es <- rnorm(k, mean = mu, sd = sqrt(tau2))
    hes <- rnorm(k, es, sqrt(2/ni)) 
    
    # Monte Carlo and Adaptive Quadrature CD-Edgington
    mc <- remaeffect(hes, se, method = "MC")
    aq <- remaeffect(hes, se, method = "GAQ")
    
    return(unname(c(mc$estimate, mc$CI, aq$estimate, aq$CI)))
  }
  
  # Grid of investigated conditions
  presimlev <- expand.grid(
    k = c(3, 5, 10, 20, 50),
    I2 = c(0.0, 0.3, 0.6, 0.9),
    i = 1:1000
  )
  
  # Parallelization
  num_cores <- parallel::detectCores() - 5
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  res_MCvsGAQ <- foreach(j = seq_len(nrow(presimlev)), 
                         .errorhandling = "remove",
                         .packages = "edgemeta",
                         .combine = rbind) %dopar% {
                           one_simulation(presimlev$k[j], presimlev$I2[j])
                         }
  stopCluster(cl)
  
  # Result
  res_MCvsGAQ <- cbind(res_MCvsGAQ, presimlev)
  colnames(res_MCvsGAQ) <- c("mcp", "mcl", "mcu", "aqp", "aql", "aqu",
                             "k", "I2", "i")
  saveRDS(res_MCvsGAQ, file = "data/MCvsGAQ.RDS")
} else {
  res_MCvsGAQ <- readRDS("data/MCvsGAQ.RDS")  
}   

# Compute bias, coverage and differences between pointe estimates and confidence
# interval limits (true average effect = -0.3)
sum_MCvsGAQ <- res_MCvsGAQ |>
  mutate(
    mc_cover = mcl <= -0.3 & mcu >= -0.3,
    aq_cover = aql <= -0.3 & aqu >= -0.3,
    mc_bias = mcp + 0.3,
    aq_bias = aqp + 0.3
  ) |> 
  group_by(k, I2) |> 
  summarise(
    mcbias        = mean(mc_bias), 
    mc_bias_se    = sd(mc_bias) / sqrt(n()),
    aqbias        = mean(aq_bias), 
    aq_bias_se    = sd(aq_bias) / sqrt(n()),
    mccover       = mean(mc_cover), 
    mc_cover_se   = sqrt(mccover * (1 - mccover) / n()),
    aqcover       = mean(aq_cover), 
    aq_cover_se   = sqrt(aqcover * (1 - aqcover) / n()),
    diffp         = mean(mcp - aqp),
    diff_p_se     = sd(mcp - aqp) / sqrt(n()),
    difflower     = mean(mcl - aql), 
    diff_lower_se = sd(mcl - aql) / sqrt(n()),
    diffupper     = mean(mcu - aqu), 
    diff_upper_se = sd(mcu - aqu) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    mcbias    = sprintf("%.3f [%.3f]", mcbias, mc_bias_se),
    aqbias    = sprintf("%.3f [%.3f]", aqbias, aq_bias_se),
    mccover   = sprintf("%.3f [%.3f]", mccover, mc_cover_se),
    aqcover   = sprintf("%.3f [%.3f]", aqcover, aq_cover_se),
    diffp     = sprintf("%.3f [%.3f]", diffp, diff_p_se),
    difflower = sprintf("%.3f [%.3f]", difflower, diff_lower_se),
    diffupper = sprintf("%.3f [%.3f]", diffupper, diff_upper_se),
  ) |>
  dplyr::select(k, I2, mcbias, aqbias, mccover, aqcover, diffp, difflower, diffupper)

# Helper functions xtables
fm <- function(dt, metric) {
  dt |>
    dplyr::select(k, I2, !!sym(metric)) |>
    pivot_wider(names_from = I2, values_from = !!sym(metric)) |>
    arrange(k) |>
    mutate(k = as.character(k)) |>
    rename_with(~ paste0(as.numeric(.x) * 100, "%"), -k)
}

header <- function(dt, section_name) {
  hr <- as.list(c(section_name, "", "", "", ""))
  i2r <- as.list(c("Studies ($k$)", "$\\iota^2$ = 0\\%", "30\\%", "60\\%", "90\\%"))
  rbind(hr, i2r, dt)
}

# Table with mean differences in point estimates and confidence interval limits
# between MC and GAQ CD-Edgington
meandiff_MCGAQ <- xtable(bind_rows(
  header(fm(sum_MCvsGAQ, "diffp"),  "Estimate (MC - GAQ)"),
  header(fm(sum_MCvsGAQ, "difflower"), "Lower 95\\% CI (MC - GAQ)"),
  header(fm(sum_MCvsGAQ, "diffupper"), "Upper 95\\% CI (MC - GAQ)")),
  align = c("l", "l", "r", "r", "r", "r"),
  caption = paste("Mean differences (with Monte Carlo standard errors) in point",
                  "estimates and 95\\% confidence interval limits between Monte",
                  "Carlo sampling and global adaptive quadrature integration for",
                  "the CD-Edgington estimator."),
  label = "tab:MCvsGAQ")

addtorow <- list()
addtorow$pos <- list(nrow(meandiff_MCGAQ))
addtorow$command <- paste(" \n\\multicolumn{5}{l}{\\footnotesize CI = confidence",
                          "interval, GAQ = global adaptive quadrature, MC =",
                          "Monte Carlo.} \\\\ \n")

print(meandiff_MCGAQ,
      include.rownames = FALSE, 
      include.colnames = FALSE,
      caption.placement = "top",
      hline.after = c(0, 2, 7, 14, 21),
      booktabs = TRUE,
      sanitize.text.function = identity,
      add.to.row = addtorow)

# Table with bias and confidence interval coverage for MC and GAQ CD-Edgington
bias_coverage_MCGAQ <- xtable(bind_rows(
  header(fm(sum_MCvsGAQ, "mcbias"),  "MC: Bias"), 
  header(fm(sum_MCvsGAQ, "aqbias"),  "GAQ: Bias"),
  header(fm(sum_MCvsGAQ, "mccover"), "MC: 95\\% CI coverage"),
  header(fm(sum_MCvsGAQ, "aqcover"), "GAQ: 95\\% CI coverage")),
  align = c("l", "l", "r", "r", "r", "r"),
  caption = paste("Bias of point estimators and coverage of 95\\% confidence",
                  "intervals (with Monte Carlo standard errors) for Monte Carlo",
                  "sampling and global adaptive quadrature integration for the",
                  "CD-Edgington estimator."),
  label = "tab:MCvsGAQ_bc")  

addtorow <- list()
addtorow$pos <- list(nrow(bias_coverage_MCGAQ))
addtorow$command <- paste(" \n\\multicolumn{5}{l}{\\footnotesize CI = confidence",
                          "interval, GAQ = global adaptive quadrature, MC =",
                          "Monte Carlo.} \\\\ \n")

print(bias_coverage_MCGAQ, 
      include.rownames = FALSE,
      include.colnames = FALSE,
      caption.placement = "top",
      hline.after = c(0, 2, 7, 14, 21, 28),
      booktabs = TRUE,
      sanitize.text.function = identity,
      add.to.row = addtorow)


## Web Appendix C: Continued Analysis of Corticosteroids and Mortality in
## Hospitalized COVID-19 Patients
## Computation
##------------------------------------------------------------------------------
es_covid <- data_covid_example$logOR  # study estimates
se_covid <- data_covid_example$logSE  # standard errors

# Random-effects meta-analysis with HTS prediction interval
meta_covid <- metagen(TE = es_covid, seTE = se_covid, method.tau = "PM", 
                      method.random.ci = "HK", method.predict = "HTS")

# Edgington (Held et al., 2025)
cm_covid <- confMeta(es_covid, se_covid, heterogeneity = "additive",
                     tau2 = meta_covid$tau2, conf_level = 0.95,
                     fun = p_edgington, 
                     fun_name = "Edgington  (one-sided input)",
                     input_p = "greater")

# CD-Edgington estimator (Monte Carlo)
CDedge <- remaeffect(es_covid, se_covid, level.ci = 0.95, seed = 982)

# Reproducibility under Monte Carlo sampling (matching manuscript)
set.seed(021098)

# PCD-Fixed (Edgington-based)
pdcovid_fixed <- PredDist(es_covid, se_covid, method = "PCD-fixed", 
                          method.tau2 = "PM")
pdcovid_fixed_pi <- pdcovid_fixed$PI

# PCD-Simplified (Edgington-based)
pdcovid_simplified <- PredDist(es_covid, se_covid, method = "PCD-simplified", 
                               method.tau2 = "PM")
pdcovid_simplified_pi <- pdcovid_simplified$PI

# PCD-Full (Edgington-based)
pdcovid_full <- PredDist(es_covid, se_covid, method = "PCD-full")
pdcovid_full_pi <- pdcovid_full$PI

# Higgins-Thompson-Spiegelhalter predictive distribution
pdhts <- rt(1e7, df = length(es_covid) - 2) * # Student's t
  sqrt(meta_covid$seTE.random^2 + meta_covid$tau2) + # scale
  meta_covid$TE.random # shift

# Reconstruct predictive distribution by Nagashima et al.
# * draw tau2 from exact confidence distribution of standard Q
pima_tau2 <- pimeta::pima(es_covid, se_covid, seed = 4833, B = 1e5)$rnd

# * Generate samples theta_new based on drawn samples of tau2
pima_sample_tn <- function(es, se, tau2) {
  
  # IVW estimator
  w_i <- 1 / (se^2 + tau2)
  mu_ivw <- sum(w_i * es) / sum(w_i)
  
  # Hartung-Knapp standard error
  q <- sum(w_i * (es - mu_ivw)^2) / (length(es) - 1)
  se_hk <- sqrt(q / sum(w_i))
  
  # Samples
  Z <- rnorm(1)
  t <- rt(1, length(es) - 1)
  
  # Theta_new = IVW + Z * tau - T * seHKSJ(IVW)
  mu_ivw + Z * sqrt(tau2) - t * se_hk
}
pdnnf <- sapply(pima_tau2, function(tau2) pima_sample_tn(es_covid, se_covid, tau2))

# Store samples from predictive distributions in list
list_pd <- list(pdhts, 
                pdnnf, 
                pdcovid_fixed$samples[, "theta_new"],
                pdcovid_simplified$samples[, "theta_new"],
                pdcovid_full$samples[, "theta_new"])

# Confidence of theta_new > 0
gr0 <- round(sapply(list_pd, function(x) mean(x >= 0)),3)

## Web Table 3: Summary table COVID-19 Data
## -----------------------------------------------------------------------------
summary_table_covid <- data_covid_example |>
  mutate(
    `OR` = sprintf("%.2f", exp(logOR)),
    `95% CI` = sprintf("%.2f to %.2f", lower, upper) 
  ) |>
  dplyr::select(Study = name, Cx = steroids, `No Cx` = nosteroids, OR, `95% CI`) |>
  xtable(
    caption = paste(
      "Summary of seven randomized controlled trials on corticosteroids and mortality",
      "in hospitalized COVID-19 patients \\citep{who2020corticosteroids}."
    ),
    label = "tab:coviddata",
    align = c("l", "l", "r", "r", "r", "r"),
    digits = c(0, 0, 0, 2, 2, 0)
  )

addtorow <- list()
addtorow$pos <- list(-1, nrow(data_covid_example))
addtorow$command <- c(
  paste("\\hline \\\\[-1em] \\multicolumn{1}{l}{} & \\multicolumn{2}{c}",
        "{Deaths / Patients} & \\multicolumn{2}{c}{} \\\\ \n"),
  paste("\\multicolumn{5}{l}{\\scriptsize CI = confidence interval,",
        "OR = odds ratio, Cx = Corticosteroids.} \\\\ \n")
)

print(summary_table_covid,
      include.rownames = FALSE,
      caption.placement = "top",
      hline.after = c(0), 
      booktabs = TRUE,
      add.to.row = addtorow)

## Web Table 4: Point estimates and confidence intervals
## -----------------------------------------------------------------------------
estimates_covid <- data.frame(
  Method = c("Hartung--Knapp--Sidik--Jonkman", "Edgington", "CD-Edgington"),
  Estimate = c(meta_covid$TE.random, cm_covid$p_max[,"x"], CDedge$estimate),
  Lower = c(meta_covid$lower.random, cm_covid$joint_cis[1], CDedge$CI[1]),
  Upper = c(meta_covid$upper.random, cm_covid$joint_cis[2], CDedge$CI[2])
  ) |>
  rowwise() |>
  mutate(
    sk = skew_pi(c(Lower, Upper), Estimate),
    sk = case_when(round(sk, 4) == 0 ~ 0,  TRUE ~ sk),
    ci = paste0(sprintf("%.2f",Lower), "  to  ", sprintf("%.2f",Upper)),
    w = Upper - Lower
  ) |>
  bind_cols(pval = biostatUZH::formatPval(c(
    meta_covid$pval.random, cm_covid$p_0[2], CDedge$pval
  ))) |>
  dplyr::select(Method, Estimate, ci, w, sk, pval) |>
  rename(
    "Method" = "Method", "95\\% CI" = "ci", "Skewness" = "sk",
    "Width" = "w", "p-value" = "pval"
  ) |> 
  xtable(caption = paste(
    "Point estimates and 95\\% confidence intervals (CI) for the average treatment",
    "effect (log odds ratio), based on seven randomized controlled trials",
    "investigating the association between corticosteroids and mortality in",
    "hospitalized COVID-19 patients \\citep{who2020corticosteroids}. Two-sided",
    "$p$-values are reported for testing $H_0: \\mu = 0$ against $H_1: \\mu \\neq 0$."
  ),
  align = c("l", "l", "r", "r", "r", "r", "r"),
  label = "tab:covidpointestimates")

print(estimates_covid,
      include.rownames = FALSE,
      caption.placement = "top",
      hline.after = c(-1,0),
      booktabs = TRUE,
      sanitize.text.function = identity)

## Web Figure 3: Monte Carlo confidence distributions of tau2 and mu
## -----------------------------------------------------------------------------
p_ctau2 <- data.frame(mc = CDedge$cd_tau2) |>
  ggplot(aes(x = mc)) + 
  geom_histogram(bins = 100, aes(y = after_stat(density)), 
                 fill = "gray95", color = "gray20") +
  geom_line(inherit.aes = FALSE, data = data.frame(
    x = seq(-0.02, 1.5, l = 100),
    y = cdtau2(seq(-0.02, 1.5, l = 100), es_covid, se_covid)),
    aes(x = x, y = y, color = "Analytic density"), linewidth = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0,7.5)) +
  scale_color_manual(
    values = c("Analytic density" = "royalblue"),
    labels = c("Analytic density\n(change of variables)"),
    name = NULL
  ) +
  labs(y = expression(c~(tau^2)), x = expression(tau^2), title = "A") +
  scale_y_continuous(expand = c(0,0.025)) +
  theme_dk() +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "gray30", 
                                     linewidth = 0.3),
    legend.margin = margin(2,2,2,2),
    legend.key.height = unit(0.4, "cm"),
    legend.text = element_text(size = 10)
  )

p_cmu <- data.frame(mc = CDedge$cd_mu) |>
  ggplot(aes(x = mc)) + 
  geom_histogram(bins = 50, aes(y = after_stat(density)), 
                 fill = "gray95", color = "gray20") +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(0, 2)) +
  geom_vline(aes(
      xintercept = CDedge$estimate,
      linetype = "Point estimate",
      color = "Point estimate"
    ),
    linewidth = 0.7
  ) +
  geom_vline(
    aes(
      xintercept = CDedge$CI[1],
      linetype = "95% CI",
      color = "95% CI"
    ),
    linewidth = 0.6
  ) +
  geom_vline(
    aes(
      xintercept = CDedge$CI[2],
      linetype = "95% CI",
      color = "95% CI"
    ),
    linewidth = 0.6
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Point estimate" = "solid",
      "95% CI" = "dashed"
    )
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Point estimate" = "#1b9e77",
      "95% CI" = "#d95f02"
    )
  ) +
  labs(y = lab_fmu, x = lab_mu, title = "B") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.0)) +
  theme_dk() +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "gray30",
                                     linewidth = 0.3),
    legend.spacing.y = unit(0.05, "cm"), 
    legend.key.height = unit(0.25, "cm"),
    legend.margin = margin(2, 2, 2, 2),
    legend.text = element_text(size = 10)
  )

p_ctau2 + p_cmu

## Web Table 5: Prediction intervals table
##------------------------------------------------------------------------------
# Function to find confidence/probability of theta_new > 0 for R meta methods
find_conf <- function(mm) {
  x <- optimize(function(lpi) {
    varest <- "PM"
    if (mm %in% c("KR", "KR-PR")) varest <- "REML"
    
    me <- metagen(es_covid, se_covid, random = T,
                  method.predict = mm,
                  method.tau = varest,
                  level.predict = lpi)
    abs(me$upper.predict)
  }, lower = 0, upper = 1)$minimum 
  
  (1 - x) / 2
}

meta_pmethods <- c("V", "HTS", "HK", "HK-PR", "KR", "KR-PR", "NNF", "S")
gr0_meta <- unlist(lapply(meta_pmethods, find_conf))
gr0_meta[2] <- gr0[1] # to match sampling (HTS)
gr0_meta[7] <- gr0[2] # to match sampling (NNF/bootstrap)

# 95% prediction intervals for all R meta methods
meta_combos <- expand.grid(method = meta_pmethods, level = 0.95)
covid_metapis <- as.data.frame(
  do.call(rbind, lapply(1:nrow(meta_combos), function(i) {
    method <- meta_combos$method[i]
    level  <- meta_combos$level[i]
    
    varest <- "PM"
    if (method %in% c("KR", "KR-PR")) varest <- "REML"
    
    m <- metagen(
      TE = es_covid,
      seTE = se_covid,
      method.tau = varest,
      B = 100000,
      method.predict = method,
      level.predict = level
    )
    
    data.frame(
      method = method,
      estimate = m$TE.random,
      pi.lower = m$lower.predict,
      pi.upper = m$upper.predict
    )
  }))
)

# Median for NNF by MC
covid_metapis[covid_metapis$method == "NNF", "estimate"] <- median(list_pd[[2]])

# Edgington's PCD distributions 
PCD_intervals <- data.frame(
  method = c("FixedTau2", "SimplifiedCD", "FullCD"),
  estimate = sapply(list_pd[3:5], median),
  pi.lower = c(pdcovid_fixed_pi[1], pdcovid_simplified_pi[1], pdcovid_full_pi[1]),
  pi.upper = c(pdcovid_fixed_pi[2], pdcovid_simplified_pi[2], pdcovid_full_pi[2])
)

# Bind prediction intervals together 
covid_pi <- covid_metapis |>
  bind_rows(PCD_intervals) |>
  rowwise() |>
  mutate(
    method = factor(method,
                    levels = rev(c("S", "HTS", "V", "HK", "HK-PR", "KR", "KR-PR",
                                   "NNF", "FixedTau2", "SimplifiedCD", "FullCD")),
                    labels = rev(c("Skipka", "HTS", "HTS-Veroniki",
                                   "HTS-HK", "HTS-HK-PR", "HTS-KR", "HTS-KR-PR",
                                   "Bootstrap", "PCD-fixed",
                                   "PCD-simplified", "PCD-full")))
  )

# Display table for 95% prediction intervals
pixtabs <- xtable(
  covid_pi |>
    bind_cols(data.frame(gr = c(gr0_meta, (gr0)[3:5])))|>
    rowwise() |>
    mutate(
      sk = skew_pi(c(pi.lower, pi.upper), estimate),
      w = ifelse(pi.upper - pi.lower > 10^6, "$>10^6$",
                 as.character(sprintf("%.2f", pi.upper - pi.lower))),
      method = factor(method,
                      levels = c("PCD-full", "PCD-simplified", "PCD-fixed",
                                 "Bootstrap", "HTS-KR-PR", "HTS-KR", "HTS-HK-PR",
                                 "HTS-HK", "HTS-Veroniki", "HTS", "Skipka")),
      sk = case_when(
        round(sk, 4) == 0 ~ 0, 
        TRUE ~ sk
      ),
      estimate = ifelse(is.na(pi.lower), "NC", as.character(sprintf("%.2f", estimate))),
      gr = ifelse(is.na(pi.lower), " ", as.character(sprintf("%.3f", gr ))),
      pi = case_when(
        is.na(pi.lower) | is.na(pi.upper) ~ "NC",
        pi.lower < -1e6 & pi.upper > 1e6 ~ "$<10^6$ to $>10^6$",
        TRUE ~ paste0(sprintf("%.2f", pi.lower), " to ", sprintf("%.2f", pi.upper))
      )
    ) |>
    dplyr::select(method, estimate, pi, w, sk, gr) |>
    rename(
      "Median " = estimate,
      "Method" = method,
      "95\\% PI" = pi,
      "Skewness" = sk,
      "Width" = w,
      "Conf($\\tn \\ge 0$)" = gr
    ) |>
    arrange(Method),
  label = "tab:serepi",
  align = c("l", "l", "r", "r", "r", "r", ">{\\raggedleft\\arraybackslash}p{3cm}"),
  digits = c(0,0,2,0,2,3,3),
  caption = paste("95\\% prediction intervals from the \\texttt{meta} package",
                  "and Edgington’s PCD distributions for seven randomized trials",
                  "on corticosteroids and COVID-19 mortality",
                  "\\citep{who2020corticosteroids}. Medians and confidence",
                  "probabilities for a future effect $\\ge 0$ are also shown.",
                  "Heterogeneity is estimated via Paule--Mandel",
                  "\\citep{paule1982consensus}, except for HTS-KR(-PR), which",
                  "require REML.")
)

addtorow <- list()
addtorow$pos <- list(nrow(pixtabs))
addtorow$command <- c(
  "\n\\multicolumn{6}{l}{\\footnotesize\\makebox[0pt][l]{%
  \\parbox[t]{0.9\\textwidth}{HK = Hartung--Knapp, KR = Kenward--Roger, NC = non-convergence,
  PI = prediction interval, \\\\ PR = Partlett--Riley.}}}\\\\\n"
)

print(pixtabs,
      caption.placement = "top",
      include.rownames = FALSE,
      hline.after = c(-1,0),
      booktabs = TRUE, 
      add.to.row = addtorow,
      sanitize.text.function = identity)

## Web Appendix D: Simulation Study Details: Example of skew-normal distribution
##------------------------------------------------------------------------------
# Parameters by moment-matching
mu = -0.3
tau2 = 0.5
alpha <- - 4
d <- alpha / sqrt(1 + alpha ^ 2)
omega <- sqrt(tau2 / (1 - 2 * (d^2) / pi))
xi <- mu - omega * d * sqrt(2 / pi)
 
pars_sn <- data.frame(expand.grid(
  x = seq(-4, 2.5, length.out = 100))) |>
  rowwise() |>
  mutate(
    dens = dsn(x, xi, omega, alpha),
    label = paste0(
      'plain("SN(") * xi==', round(xi,2),
      '*","~omega==', round(omega, 2),
      '*","~alpha==', alpha,
      '*")"'
    )
  )

ggplot(pars_sn, aes(x = x, y = dens, color = alpha, group = alpha)) +
  facet_grid(~label,
             labeller = label_parsed) +
  geom_line(color = "black", linewidth = 0.5) +
  labs(y = "Density", x = "Y") +
  scale_y_continuous(expand = c(0,0.01), limits = c(0,0.65)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_dk() +
  theme(
    strip.background = element_rect(fill = "gray99", color = "black"),
    legend.position = "none"
  )


## Web Appendix E: Additional Simulation Results
## Load and pre-process results
## -----------------------------------------------------------------------------
simresi <- list.files("../results", full.names = TRUE)
dt_v <- do.call(rbind, lapply(simresi, function(i) {
  x <- readRDS(i)
})) |>
  mutate(
    k = factor(k, levels = c(3,5,10,20,50)),
    I2 = factor(I2, levels = c(0, 0.3, 0.6, 0.9),
                labels = c(0, 30, 60, 90))
  ) |>
  rename("sqe.held.a" = "seq.held.a")

# Listwise deletion
dt <- dt_v[complete.cases(dt_v), ]

## E.2 Point Estimation: Maximum bias of point estimators
## -----------------------------------------------------------------------------
# Normal effects
bmax <- min(dt |>
              filter(dist == "N") |>
              group_by(I2, k, k_large) |>
              summarise(
                bmax = mean(bias.ivw.random), .groups = "drop"
              ) |>
              pull(bmax))

# Skew normal effects
bmaxlsn <- max(dt |>
                 filter(dist == "LSN") |>
                 group_by(I2, k, k_large) |>
                 summarise(
                   bmax = mean(bias.held.a), .groups = "drop"
                 ) |>
                 pull(bmax))

## Web Table 6: Non-convergences parametric bootstrap prediction interval
## -----------------------------------------------------------------------------
dt_v[!complete.cases(dt_v), ] |>
  mutate(
    k_large = as.character(k_large),
    dist = case_when(
      dist == "N" ~ "Normal", 
      dist == "LSN" ~ "Skew-normal"
    )) |>
  group_by(k, I2, k_large, dist) |>
  summarise("Non-convergences" = n()) |>
  rename(
    "Studies" = "k",
    "Large studies" = "k_large",
    "$\\iota^2$ (\\%)" = "I2",
    "Distribution" = "dist"
  ) |> 
  xtable(label = "tab:ncnnf",
         align = c("r", "r", "r", "r", "l", "r"),
         caption = paste("Scenarios of the simulation study under which the",
                         "95\\% parametric bootstrap prediction interval did not",
                         "converge at least once.")
  )|>
  print(include.rownames = FALSE,
        sanitize.colnames.function = identity,
        caption.placement = "top",
        hline.after = c(-1,0),
        booktabs = TRUE)

## Web Figure 5: Coverage of 95% Prediction Intervals under skew-normal effects
##------------------------------------------------------------------------------
hlines <- expand.grid(
  k_large = paste0("Large~studies==", 0:2),
  I2 = "iota^2==0*'%'"
)
hlines2 <- expand.grid(
  k_large = paste0("Large~studies==", 0:2),
  I2 = c("iota^2==30*'%'", "iota^2==60*'%'", "iota^2==90*'%'")
)

simplot(dt, "LSN", "pi.cvr.", nam = 1) +
  geom_hline( 
    data = hlines,
    aes(yintercept = 1.00),
    linetype = "dashed", size = 0.25
  ) + 
  geom_hline(
    data = hlines2,
    aes(yintercept = 0.95),
    linetype = "dashed", size = 0.25
  ) +
  labs(y = "Coverage of 95% prediction intervals \u00B1 MCSE") +
  scale_y_continuous(labels = scales::label_percent(),
                     breaks = seq(0.8, 1, by = 0.05),
                     expand = c(0.02, 0))


## Web Figures 6 to 12: 95% PI coverage histograms
## -----------------------------------------------------------------------------
pcvr(0, NULL) + labs(y = "Density", x = "Coverage probability") 
pcvr(30, "N")         
pcvr(30, "LSN")         
pcvr(60, "N")         
pcvr(60, "LSN")         
pcvr(90, "N")          
pcvr(90, "LSN")          

## Web Figure 13: Skewness agreement of PIs with effect estimates
## -----------------------------------------------------------------------------
# Skew-normal effects
plotcor(dt = dt, distr = "LSN", stw_var = "pi.sk.", var = sk.hes, nam = 1) +
  labs(y = "Pearson correlation \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## Web Figures 14 and 15: Skewness sign agreement of PIs with effect estimates 
## (Cohens Kappa)
## -----------------------------------------------------------------------------
# Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.hes, nam = 1) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

# Skew-Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.hes, nam = 1) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## Web Figures 16 and 17: Width of 95% Prediction Intervals
##------------------------------------------------------------------------------
# Normal effects
simplot(dt = dt, distr = "N", stw_var = "pi.w.", nam = 1) +
  labs(y = "Width of 95% prediction intervals \u00B1 MCSE") +
  scale_y_log10()

# Skew-Normal effects
simplot(dt = dt, distr = "LSN", stw_var = "pi.w.", nam = 1) +
  labs(y = "Width of 95% prediction intervals \u00B1 MCSE") +
  scale_y_log10()

## Web Figures 18 and 19: CRPS
##------------------------------------------------------------------------------
# Normal effects
simplot(dt = dt, distr = "N", stw_var = "crps.", nam = 2) +
  labs(y = "CRPS \u00B1 MCSE")

# Skew-Normal effects
simplot(dt = dt, distr = "LSN", stw_var = "crps.", nam = 2) +
  labs(y = "CRPS \u00B1 MCSE")

## Web Figures 20 and 21: Computation time of 95% prediction intervals
##------------------------------------------------------------------------------
# Normal effects
simplot(dt, "N", "t.", 5) +
  labs(y = "Computation time (seconds) \u00B1 MCSE") +
  scale_y_continuous(trans = "log2") 

# Skew-normal effects
simplot(dt, "LSN", "t.", 5) +
  labs(y = "Computation time (seconds) \u00B1 MCSE") +
  scale_y_continuous(trans = "log2") 

## Web Figures 22 and 23: Coverage 95% confidence intervals
##------------------------------------------------------------------------------
# Normal effects 
simplot(dt, "N", "ci.cvr.", nam = 4) +
  labs(y = "Coverage of 95% confidence intervals \u00B1 MCSE") + 
  geom_hline(yintercept = 0.95, linetype = "dashed", size = 0.25) +
  scale_y_continuous(labels = scales::label_percent(),
                     breaks = c(0.7, 0.8, 0.9, 1))

# Skew-normal effects
simplot(dt, "LSN", "ci.cvr.", nam = 4) +
  labs(y = "Coverage of 95% confidence intervals \u00B1 MCSE") +
  geom_hline(yintercept = 0.95, linetype = "dashed", size = 0.25) +
  scale_y_continuous(labels = scales::label_percent())

## Web Figures 24 and 25: Width of 95% confidence intervals
## -----------------------------------------------------------------------------
# Normal effects
simplot(dt = dt, distr = "N", stw_var = "ci.w.", nam = 4) +
  labs(y = "Width of 95% confidence intervals \u00B1 MCSE")

# Skew-normal effects
simplot(dt = dt, distr = "LSN", stw_var = "ci.w.", nam = 4) +
  labs(y = "Width of 95% confidence intervals \u00B1 MCSE")

## Web Figures 26 and 27: Skewness agreement of CIs with effect estimates
## -----------------------------------------------------------------------------
# Normal effects
plotcor(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.hes, nam = 2) +
  labs(y = "Pearson correlation \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

# Skew-normal effects
plotcor(dt = dt, distr = "LSN", stw_var = "ci.sk.", var = sk.hes, nam = 2) +
  labs(y = "Pearson correlation \u00B1 MCSE")  +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## Web Figures 28 and 29: Skewness sign agreement of CIs with effect estimates
## (Cohens kappa)
## -----------------------------------------------------------------------------
# Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.hes, nam = 2) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

# Skew-Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.hes, nam = 2) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## Web Figures 30 and 31: Bias point estimator
##------------------------------------------------------------------------------
# Normal effects
simplot(dt = dt, distr = "N", stw_var = "bias.", nam = 3) +
  labs(y = "Bias \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  scale_y_continuous(limits = c(-0.05, 0.05))

# Skew-normal effects
simplot(dt = dt, distr = "LSN", stw_var = "bias.", nam = 3) +
  labs(y = "Bias \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## Web Figures 32 and 33: MSE point estimator
##------------------------------------------------------------------------------
# Normal effects
simplot(dt, "N", "sqe.", nam = 3) +
  labs(y = "MSE \u00B1 MCSE")+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

# Skew-normal effects
simplot(dt, "LSN", "sqe.", nam = 3) +
  labs(y = "MSE \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## Web Table 8: Computation time of Monte Carlo CD-Edgington
##------------------------------------------------------------------------------
if (!file.exists("data/timeCDEdgington.RDS")) {
  
  one_simulation <- function(k) {

    ni <- 50
    mu <- -0.3
    I2 <- 0.3
    
    tau2 <- 1 / k * sum(2 / ni) * (I2 / (1 - I2))
    se <- sqrt(rchisq(k, df = 2 * (ni - 1)) * ((ni - 1) * ni)^(-1))
    
    es <- rnorm(k, mean = mu, sd = sqrt(tau2))
    hes <- rnorm(k, es, sqrt(2/ni)) 
    
    return(system.time(remaeffect(hes, se, method = "MC"))["elapsed"])
  }
  
  # Grid of investigated conditions
  fl <- expand.grid(k = c(3, 5, 10, 20, 50), i = 1:1000)
  
  # Parallelization
  num_cores <- parallel::detectCores() - 5
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  res <- foreach(j = seq_len(nrow(fl)), 
                 .errorhandling = "remove",
                 .packages = "edgemeta",
                 .combine = rbind) %dopar% {
                   one_simulation(fl$k[j])
                 }
  stopCluster(cl)
  saveRDS(cbind(fl, res), file = "data/timeCDEdgington.RDS")
} 

# Summarize computation time in table 
readRDS("data/timeCDEdgington.RDS") |>
  group_by(k) |>
  summarise(
    t = paste0(sprintf("%.2f", mean(elapsed)),
               " [",
               sprintf("%.2f", sd(elapsed) / sqrt(n())),
               "]")) |>
  mutate(k = factor(k)) |>
  rename("Studies ($k$)" = "k", "Runtime (s) [MCSE]" = "t") |>
  xtable(
    caption = paste("Mean computation time (seconds; s) with Monte Carlo",
                    "standard error (MCSE) of the Monte Carlo CD-Edgington",
                    "estimator across varying numbers of studies, based on 1000",
                    "simulation iterations."),
    label = "tab:timeCDEdgington",
    align = c("l", "l", "r")
  ) |>
  print(include.rownames = FALSE,
        caption.placement = "top",
        hline.after = c(-1, 0), 
        sanitize.text.function = identity,
        booktabs = TRUE)

## Web Appendix F: Weighted Edgington Methods
## Computation
## -----------------------------------------------------------------------------
# Weights
w_unweighted <- rep(1, length(es_covid))
w_se <- 1 / se_covid
w_se2 <- 1 / se_covid^2

# 95% confidence intervals
est_unweighted <- CDedge
est_w_se <- remaeffect(es_covid, se_covid, method = "MC", w = w_se)
est_w_se2 <- remaeffect(es_covid, se_covid, method = "MC", w = w_se2)

# 95% prediction intervals
pi_unweighted <- pdcovid_full
pi_w_se <- PredDist(es_covid, se_covid, method = "PCD-full", w = w_se)
pi_w_se2 <- PredDist(es_covid, se_covid, method = "PCD-full", w = w_se2)

# Display prediction and confidence intervals
# ------------------------------------------------------------------------------
n_pi  <- 5 # number of prediction methods
n_ci  <- 4 # number of estimation methods
gap1  <- 0.6  # vertical gap between study estimates and estimation methods
gap2  <- 0.9  # vertical gap between estimation and prediction methods
study_gap <- 0.55 # smaller gaps between studies
xmin_plot <- -2.5 # lower x-axis-limit
xmax_plot <- 2 # upper x-axis limit

# Add number of patients to data
data_covid_example <- data_covid_example |>
  mutate(
    lab_N = as.character(
        as.integer(sub(".*/", "", steroids)) + as.integer(sub(".*/", "", nosteroids))
      )
    )

# Helper for preparing plotting data
mk <- function(label, estimate = NA_real_, lower, upper) {
  tibble(label, method = label, estimate, lower, upper)
}

# KDE helper
kde <- function(label, x) { 
  d <- density(x)
  tibble(label, method = label, x = d$x, d = d$y) 
}

# Labels on y-axis
label_order <- c(
  "PCD-full (w = 1 / SE2)", 
  "PCD-full (w = 1 / SE)", 
  "PCD-full (unweighted)",
  "Parametric bootstrap", "HTS",
  "CD-Edgington (w = 1 / SE2)", "CD-Edgington (w = 1 / SE)", 
  "CD-Edgington (unweighted)",
  "Hartung-Knapp-\nSidik-Jonkman",
  rev(data_covid_example$name)
)
label_expr <- as.expression(c(
  expression("PCD-full (w=" * 1 / hat(sigma)[i]^2 * ")"),
  expression("PCD-full (w=" * 1 / hat(sigma)[i] * ")"),
  "PCD-full (unweighted)",
  "Parametric bootstrap",
  "Higgins-Thompson-\nSpiegelhalter",
  expression("CD-Edgington (w=" * 1 / hat(sigma)[i]^2 * ")"),
  expression("CD-Edgington (w=" * 1 / hat(sigma)[i] * ")"),
  "CD-Edgington (unweighted)",
  "Hartung-Knapp-\nSidik-Jonkman",
  rev(data_covid_example$name)
))

# Combine study estimates, estimation and prediction results
dt_weighted_example <- bind_rows(
  data_covid_example |> transmute(label = name, method = "Study", estimate = logOR,
                          lower = logOR - 1.96 * logSE, upper = logOR + 1.96 * logSE,
                          lab_N = lab_N),
  mk("Hartung-Knapp-\nSidik-Jonkman",
     meta_covid$TE.random,  meta_covid$lower.random, meta_covid$upper.random),
  
  mk("CD-Edgington (unweighted)",
     est_unweighted$estimate, est_unweighted$CI[1], est_unweighted$CI[2]),
  
  mk("CD-Edgington (w = 1 / SE)",
     est_w_se$estimate, est_w_se$CI[1], est_w_se$CI[2]),
  
  mk("CD-Edgington (w = 1 / SE2)",  
     est_w_se2$estimate, est_w_se2$CI[1], est_w_se2$CI[2]),
  
  mk("HTS",
     lower = meta_covid$lower.predict, upper = meta_covid$upper.predict),
  
  mk("Parametric bootstrap",
     lower = as.numeric(covid_pi[covid_pi$method == "Bootstrap", "pi.lower"]),
     upper = as.numeric(covid_pi[covid_pi$method == "Bootstrap", "pi.upper"])),
  
  mk("PCD-full (unweighted)", 
     lower = pi_unweighted$PI[1], upper = pi_unweighted$PI[2]),
  
  mk("PCD-full (w = 1 / SE)", 
     lower = pi_w_se$PI[1], upper = pi_w_se$PI[2]),
  
  mk("PCD-full (w = 1 / SE2)", 
     lower = pi_w_se2$PI[1], upper = pi_w_se2$PI[2])
) |>
  mutate(label = factor(label, label_order), y = as.numeric(label),
         upper_plot = if_else(method == "Study", pmin(upper, xmax_plot), upper),
         lower_plot = if_else(method == "Study", pmax(lower, xmin_plot), lower),
         clip_right = method == "Study" & upper > xmax_plot,
         y_plot = case_when(
           y <= n_pi ~ y,
           y <= n_pi + n_ci ~ y + gap1,
           TRUE ~
             n_pi + n_ci + gap1 + gap2 +
             (y - (n_pi + n_ci)) * study_gap
         ))

# Kernel density estimates of confidence and predictive distributions
dens_weighted_example <- bind_rows(
  kde("CD-Edgington (unweighted)", est_unweighted$cd_mu),
  kde("CD-Edgington (w = 1 / SE)", est_w_se$cd_mu),
  kde("CD-Edgington (w = 1 / SE2)", est_w_se2$cd_mu),
  tibble(label = "Hartung-Knapp-\nSidik-Jonkman", 
         method = "Hartung-Knapp-\nSidik-Jonkman",
         x = seq(meta_covid$lower.random - 1,
                 meta_covid$upper.random + 1,
                 length.out = 512),
         d = dt((x - meta_covid$TE.random) / meta_covid$seTE.random,
                df = length(es_covid) - 1) / meta_covid$seTE.random),
  kde("PCD-full (unweighted)", pi_unweighted$samples[, "theta_new"]),
  kde("PCD-full (w = 1 / SE)", pi_w_se$samples[, "theta_new"]),
  kde("PCD-full (w = 1 / SE2)", pi_w_se2$samples[, "theta_new"]),
  kde("HTS", pdhts),
  kde("Parametric bootstrap", pdnnf)
) |>
  mutate(label = factor(label, label_order), y = as.numeric(label)) |>
  group_by(label) |>
  mutate(d = ifelse(d > 0.02, d / max(d) * 0.65 + 0.01, NA)) |>
  ungroup() |>
  mutate(
    y_plot = case_when(
      y <= n_pi ~ y,
      y <= n_pi + n_ci ~ y + gap1,
      TRUE ~
        n_pi + n_ci + gap1 + gap2 +
        (y - (n_pi + n_ci)) * study_gap
    )
  )

y_breaks <- rev(dt_weighted_example$y_plot)
y_lims <- range(y_breaks) + c(-0.2, 0.2)

# Forest plot
p_forest <- ggplot(dt_weighted_example, aes(y = y_plot, color = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = .4, color = "grey45") +
  geom_line(data = dens_weighted_example, 
            aes(x, y_plot + d, group = label, color = method),
            inherit.aes = FALSE, linewidth = .75) +
  geom_errorbarh(aes(xmin = lower_plot, xmax = upper_plot),
                 height = .2, linewidth = .5) +
  geom_segment(data = subset(dt_weighted_example, clip_right),
               aes(x = xmax_plot - 0.2, xend = xmax_plot-0.05,
                   y = y_plot, yend = y_plot),
               inherit.aes = FALSE, arrow = arrow(length = unit(.13, "cm")),
               color = "black", linewidth = .5) +
  geom_point(aes(x = estimate), size = 1.5, na.rm = TRUE) +
  scale_y_continuous(breaks = y_breaks, labels = label_expr, limits = y_lims,
                     expand = c(0,0)) +
  scale_x_continuous(breaks = seq(xmin_plot, xmax_plot, 1), expand = c(0,0)) +
  coord_cartesian(xlim = c(xmin_plot, xmax_plot-0.05)) +
  scale_color_manual(values = c(
    "Study" = "black",
    "CD-Edgington (unweighted)" = "#6C91BF",
    "CD-Edgington (w = 1 / SE)" = "#88C0A9",
    "CD-Edgington (w = 1 / SE2)" = "#ACD2D9",
    "Hartung-Knapp-\nSidik-Jonkman" = "#F1C06E",
    "PCD-full (unweighted)" = "#6C91BF",   
    "PCD-full (w = 1 / SE)" = "#88C0A9",   
    "PCD-full (w = 1 / SE2)" = "#ACD2D9",  
    "HTS" = "#F1C06E",
    "Parametric bootstrap" = "#B56576"
  ), guide = "none") +
  labs(x = "Log odds ratio", y = NULL) +
  theme_dk() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey95"),
    panel.border = element_rect(colour = "white", fill = NA, linetype = 1),
    axis.line = element_line(color = "black"),         
    axis.line.x.top = element_blank(),                 
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  ) 

# Table with estimates and confidence intervals
p_table <- dt_weighted_example |>
  mutate(
    lab_est = ifelse(is.na(estimate), "", sprintf("%.2f", estimate)),
    lab_int = ifelse(is.na(lower) | is.na(upper), "",
                     sprintf("[%.2f, %.2f]", lower, upper))
  ) |>
  ggplot(aes(y = y_plot)) +
  geom_text(aes(x = 1.00, label = lab_est), hjust = 1, size = 3) +
  geom_text(aes(x = 2.85, label = lab_int), hjust = 1, size = 3) +
  geom_text(aes(x = 4.00, label = lab_N), hjust = 1, size = 3) +
  scale_x_continuous(
    limits = c(0.35, 4.35),
    breaks = c(0.9, 2.25, 3.9),
    labels = c("logOR", "95% CI / PI", "N"),
    position = "top"
  ) +
  scale_y_continuous(
    breaks = y_breaks,
    limits = y_lims,
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = y_lims, clip = "off") +
  theme_void() +
  theme(
    axis.text.x.top = element_text(size = 10, face = "bold", margin = margin(b = 6)),
    plot.margin = margin(5.5, 5.5, 5.5, 8)
  )

# Align
p_forest | p_table

# The following warning message
# <<Removed 3229 rows containing missing values or values outside the scale
# range (`geom_line()`)>>
# occurs since densities are not plotted on the entire x-range, but where greater
# than some treshold (see line 1373)

# The following warning message 
# << Removed 9 rows containing missing values or values outside the scale
# range (`geom_text()`). >>
# occurs since N is not displayed for the estimation and prediction methods

# Display 95% confidence intervals and prediction interval tables
# ------------------------------------------------------------------------------
dt_weighted_tab <- dt_weighted_example |>
  mutate(
    estimate = case_match(
      label,
      "HTS" ~ covid_metapis[covid_metapis$method == "HTS", "estimate"],
      "Parametric bootstrap" ~ covid_metapis[covid_metapis$method == "NNF", "estimate"],
      "PCD-full (unweighted)" ~ median(pi_unweighted$samples[, "theta_new"]),
      "PCD-full (w = 1 / SE)" ~ median(pi_w_se$samples[, "theta_new"]),
      "PCD-full (w = 1 / SE2)" ~ median(pi_w_se2$samples[, "theta_new"]),
      .default = estimate
    )
  ) |>
  filter(method != "Study") |>
  rowwise() |>
  mutate(
    sk = skew_pi(c(lower, upper), estimate),
    sk = sprintf("%.3f", case_when(
        round(sk, 4) == 0 ~ 0, 
        TRUE ~ sk)),
    interval = paste0(sprintf("%.2f",lower), "  to  ", sprintf("%.2f", upper)),
    w = upper - lower
  ) |>
  bind_cols(pval = biostatUZH::formatPval(c(
     meta_covid$pval.random, est_unweighted$pval, est_w_se$pval, est_w_se2$pval, rep(NA, 5)
    ))) |>
  dplyr::select(label, estimate, interval, w, sk, pval) |>
  rename(
    "Method" = "label",
    "Skewness" = "sk",
    "Width" = "w",
    "p-value" = "pval"
  ) 

# Confidence-interval table
print(
  xtable(
    dt_weighted_tab |> 
      rename("Estimate" = "estimate",
             "95\\% CI" = "interval") |>
      filter(Method %in% c(
        "Hartung-Knapp-\nSidik-Jonkman",
        "CD-Edgington (unweighted)",
        "CD-Edgington (w = 1 / SE)",
        "CD-Edgington (w = 1 / SE2)"
      )) |> 
      mutate( 
        Method = case_match(
          Method,
          "Hartung-Knapp-\nSidik-Jonkman" ~ "Hartung--Knapp--Sidik--Jonkman",
          "CD-Edgington (unweighted)" ~ "CD-Edgington (unweighted)",
          "CD-Edgington (w = 1 / SE)" ~ "CD-Edgington ($w_i = 1/\\hat\\sigma_i$)",
          "CD-Edgington (w = 1 / SE2)" ~ "CD-Edgington ($w_i = 1/\\hat\\sigma_i^2$)",
        )
      ),
    align = c("l", "l", "r", "r", "r", "r", "r"),
    caption = paste(
      "Point estimates and 95\\% confidence intervals (CI) for the average treatment",
      "effect (log odds ratio) for the corticosteroids and COVID-19 mortality",
      "meta-analysis \\citep{who2020corticosteroids}. For",
      "CD-Edgington, results are reported for unweighted, inverse-standard-error",
      "($w_i = 1/\\hat\\sigma_i$), and inverse-variance",
      "($w_i = 1/\\hat\\sigma_i^2$) weighting schemes. Results from",
      "the Hartung--Knapp--Sidik--Jonkman method are shown for reference. 
      Two-sided $p$-values are",
      "reported for testing $H_0: \\mu = 0$ against $H_1: \\mu \\neq 0$."
    ),
    label = "tab:weighted-ci-example"
  ),
  include.rownames = FALSE,
  caption.placement = "top",
  hline.after = c(-1,0), 
  booktabs = TRUE,
  sanitize.text.function = identity
)

# Prediction-interval table
conf_pi_w_se  <- suppressMessages(suppressWarnings(
  as.numeric(capture.output(conf_w_se <- conf(pi_w_se, 0))[2])
))
conf_pi_w_se2 <- suppressMessages(suppressWarnings(
  as.numeric(capture.output(conf_w_se2 <- conf(pi_w_se2, 0))[2])
))

print(
  xtable(
    dt_weighted_tab |>
      rename("Median" = "estimate",
             "95\\% PI" = "interval") |>
      filter(Method %in% c(
        "HTS",
        "Parametric bootstrap",
        "PCD-full (unweighted)",
        "PCD-full (w = 1 / SE)",
        "PCD-full (w = 1 / SE2)"
      )) |>
      dplyr::select(-`p-value`) |>
      mutate(
        "Conf($\\theta_{new} \\ge 0$)" = sprintf("%.3f", case_match(
          Method,
          "HTS" ~ gr0[1],
          "Parametric bootstrap" ~ gr0[2],
          "PCD-full (unweighted)" ~ gr0[5],
          "PCD-full (w = 1 / SE)" ~ conf_pi_w_se,
          "PCD-full (w = 1 / SE2)" ~ conf_pi_w_se2
        )),
        Method = case_match(
          Method,
          "HTS" ~
            "Higgins--Thompson--Spiegelhalter",
          "Parametric bootstrap" ~
            "Parametric bootstrap",
          "PCD-full (unweighted)" ~
            "PCD-full (unweighted)",
          "PCD-full (w = 1 / SE)" ~
            "PCD-full ($w_i = 1/\\hat\\sigma_i$)",
          "PCD-full (w = 1 / SE2)" ~
            "PCD-full ($w_i = 1/\\hat\\sigma_i^2$)",
          .default = Method
        )
      ),
    align = c("l", "l", "r", "r", "r", "r", "r"),
    caption = paste(
      "95\\% prediction intervals (PI) and medians of predictive distributions",
      "for the treatment",
      "effect in a future study (log odds ratio) for the corticosteroids and",
      "COVID-19 mortality meta-analysis \\citep{who2020corticosteroids}.",
      "PCD-full results are shown for unweighted, inverse-standard-error",
      "($w_i = 1/\\hat\\sigma_i$), and inverse-variance",
      "($w_i = 1/\\hat\\sigma_i^2$) weighting schemes.",
      "The Higgins--Thompson--Spiegelhalter and",
      "parametric bootstrap methods are included for comparison."
    ),
    label = "tab:weighted-pi-example"
  ),
  include.rownames = FALSE,
  caption.placement = "top",
  hline.after = c(-1, 0),
  booktabs = TRUE,
  sanitize.text.function = identity
)

