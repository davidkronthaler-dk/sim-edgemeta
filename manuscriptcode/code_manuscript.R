## Clear environment
##------------------------------------------------------------------------------
rm(list = ls())

## Library Packages
## -----------------------------------------------------------------------------
# the following packages are available from CRAN 
# (install with install.packages("PACKAGE"))
library(meta)
library(dplyr)
library(ggplot2)
library(patchwork)          
library(ggthemes) 
library(xtable)  
library(latex2exp) 
library(tidyr)  
library(confMeta) 
library(metafor)
library(coda) 
library(EnvStats)
library(pimeta)
library(ggnewscale)
library(scales)

# the following package is available from Github
# (install with remotes::install_github("davidkronthaler-dk/edgemeta))
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

## COVID-19 example: Random-effects meta-analysis
##------------------------------------------------------------------------------
load("data/Covid_example.Rdata")

# Extract study estimates and standard errors
es_covid <- data_covid_example$logOR  
se_covid <- data_covid_example$logSE   

# REMA with HKSJ confidence and HTS prediction interval
covid_meta <- metagen(TE = es_covid, seTE = se_covid, method.tau = "PM", 
                      method.random.ci = "HK", method.predict = "HTS")

# Edgington (Held et al., 2025)
cm_covid <- confMeta(es_covid, se_covid, heterogeneity = "additive",
                     tau2 = covid_meta$tau2, conf_level = 0.95,
                     fun = p_edgington, 
                     fun_name = "Edgington  (one-sided input)",
                     input_p = "greater")

# CD-Edgington
cd_edgington <- remaeffect(es_covid, se_covid, level.ci = 0.95, seed = 982)

# Skewness of CD-Edgington's confidence distribution
skECD <- EnvStats::skewness(cd_edgington$cd_mu)

# Confidence of average effect < 0 according to CD-Edgington
confcov <- mean(cd_edgington$cd_mu < 0)


## COVID-19 example: Forest plot (Figure 1)
## -----------------------------------------------------------------------------
adder  <- 0.4  # defines where confidence density is plotted: 95%CI +- adder
lseq   <- 400  # Length grid sequence
scaler <- 0.5  # scales confidence density height in forest plot

# CD-Edgington confidence density 
fcd_CDE <- remaeffect(es_covid, se_covid, "GAQ")$fcd
seq_CDE <- seq(cd_edgington$CI[1] - adder, cd_edgington$CI[2] + adder, l = lseq) 
dat_CDE <- data.frame(mu = seq_CDE, cd = fcd_CDE(seq_CDE)) 
dat_CDE$scaled <- dat_CDE$cd / max(dat_CDE$cd) * scaler 

# HKSJ confidence density
seq_HKSJ <- seq(covid_meta$lower.random - adder, 
                covid_meta$upper.random + adder, 
                l = lseq)
seTE_HKSJ <- covid_meta$seTE.random
TE_HKSJ   <- covid_meta$TE.random
dat_HKSJ  <- data.frame(
  mu = seq_HKSJ,
  cd = dt((seq_HKSJ - TE_HKSJ) / seTE_HKSJ, df = length(es_covid) - 1) / seTE_HKSJ
)
dat_HKSJ$scaled <- dat_HKSJ$cd / max(dat_HKSJ$cd) * scaler

# Edgington (Held et al., 2025) confidence density
fcd_EDGE <- function(mu) {
  edgemeta:::CD_cpp(mu, es_covid, sqrt(se_covid^2 + covid_meta$tau2), 
                    rep(1, length(es_covid)))
}
seq_EDGE <- seq(cm_covid$joint_cis[1] - adder, cm_covid$joint_cis[2] + adder, l = lseq)
dat_EDGE <- data.frame(mu = seq_EDGE, cd = fcd_EDGE(seq_EDGE))
dat_EDGE$scaled <- dat_EDGE$cd / max(dat_EDGE$cd) * scaler

# Data for forest plot 
dat_forest <- data_covid_example |>
  mutate(
    l = logOR - 1.96 * logSE,
    u = logOR + 1.96 * logSE,
    t = "study"
  ) |>
  bind_rows(
    data.frame(
      name = c("Edgington", "CD-Edgington", "Hartung-Knapp-\nSidik-Jonkman"),
      logOR = c(cm_covid$p_max[,"x"], cd_edgington$estimate, covid_meta$TE.random),
      l = c(cm_covid$joint_cis[1], cd_edgington$CI[1], covid_meta$lower.random),
      u = c(cm_covid$joint_cis[2], cd_edgington$CI[2], covid_meta$upper.random),
      t = c("e0", "e1", "e2") 
    )
  ) |>
  mutate(
    name = factor(name, levels = c("Hartung-Knapp-\nSidik-Jonkman", 
                                   "CD-Edgington", "Edgington", 
                                   rev(data_covid_example$name))),
    ynum = rev(c(1.0, 1.75, 2.5, seq(3.25, by = 0.4, length.out = n() - 3)))
  )

# Clip wide confidence intervals
xmax_plot <- 2
dat_forest <- dat_forest |>
  mutate(
    u_plot = pmin(u, xmax_plot),
    clipped_right = t == "study" & u > xmax_plot
  )

# Table right of forest plot (logOR, 95% CI, N)
p_table <- dat_forest |>
  mutate(
    lab_logOR = ifelse(is.na(logOR), "", sprintf("%.2f", logOR)),
    lab_logCI = ifelse(is.na(l) | is.na(u), "", sprintf("[%.2f, %.2f]", l, u)),
    lab_N = case_when(
      t == "study" ~ as.character(
        as.integer(sub(".*/", "", steroids)) + as.integer(sub(".*/", "", nosteroids))
      ),
      TRUE ~ ""
    )
  ) |>
  ggplot(aes(y = ynum)) +
  geom_text(aes(x = 1.25, label = lab_logOR), hjust = 1, size = 3) +
  geom_text(aes(x = 2.5, label = lab_logCI), hjust = 1, size = 3) +
  geom_text(aes(x = 3.15, label = lab_N), hjust = 1, size = 3) +
  scale_x_continuous(
    limits = c(0.8, 3.6),
    breaks = c(1., 2.15, 3.075),
    labels = c("logOR", "95% CI", "N"),
    position = "top"
  ) +
  theme_void() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.x.top = element_text(size = 10, face = "bold"), # optional styling
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  )

# Forest plot
p_forest <- ggplot(dat_forest, aes(x = logOR, y = ynum, color = t)) +
  geom_errorbar(aes(xmin = l, xmax = u), width = 0.05) +
  geom_segment(
    data = subset(dat_forest, clipped_right),
    aes(x = 1.95, xend = 0.99 * xmax_plot, y = ynum, yend = ynum),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.15, "cm")),
    color = "black"
  ) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_line(
    data = dat_CDE, 
    aes(x = mu, y = 1.75 + scaled),
    color = "#F1C06E", 
    linewidth = 0.8
  ) +
  geom_line(
    data = dat_HKSJ,
    aes(x = mu, y = 1 + scaled),
    color = "#90A4AE",
    linewidth = 0.8
  ) +
  geom_line(
    data = dat_EDGE,
    aes(x = mu, y = 2.5 + scaled),
    color = "#6C91DE",
    linewidth = 0.8
  ) +
  scale_color_manual(values = c("#6C91DE", "#F1C06E", "#90A4AE", "black")) +
  scale_y_continuous(breaks = dat_forest$ynum, labels = dat_forest$name) +
  scale_x_continuous(
    expand = c(0.0, 0),
    breaks = seq(-2, 2, 1)
  ) +
  coord_cartesian(xlim = c(-2, xmax_plot)) +
  labs(y = "", x = "Log odds ratio") +
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

# Align forest plot and tabular info
p_forest | p_table

## Covid-19 example: Computation of predictive distributions
##------------------------------------------------------------------------------
# Reproducibility
set.seed(021098)

# PCD-Fixed Tau2 (Edgington-based)
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

# Higgins-Thompson-Spiegelhalter (HTS)
pdhts <- rt(1e7, df = length(es_covid) - 2) * # Student's t
  sqrt(covid_meta$seTE.random^2 + covid_meta$tau2) + # scale
  covid_meta$TE.random # shift

# Reconstruct predictive distribution by Nagashima et al.
# * draw tau2 from exact confidence distribution of standard Q
pima_tau2 <- pimeta::pima(es_covid, se_covid, seed = 4833, B = 1e5)$rnd

# * Generate samples theta_new based on drawn samples of tau2
pima_sample_tn <- function(es, se, tau2) {
  
  # IVW estimator
  w_i <- 1 / (se^2 + tau2)
  mu_ivw <- sum(w_i * es) / sum(w_i)
  
  # Hartung-Knapp-Sidik-Jonkman standard error
  q <- sum(w_i * (es - mu_ivw)^2) / (length(es) - 1)
  se_hksj <- sqrt(q / sum(w_i))
  
  # Samples
  Z <- rnorm(1)
  t <- rt(1, length(es) - 1)
  
  # Theta_new = IVW + Z * tau - T * seHKSJ(IVW)
  mu_ivw + Z * sqrt(tau2) - t * se_hksj
}
pdnnf <- sapply(pima_tau2, function(tau2) pima_sample_tn(es_covid, se_covid, tau2))

# Store samples from predictive distributions in list
list_pd <- list(pdhts, pdnnf, pdcovid_fixed$samples[, "theta_new"],
                pdcovid_simplified$samples[, "theta_new"],
                pdcovid_full$samples[, "theta_new"])

# Confidence of theta_new > 0
gr0 <- round(sapply(list_pd, function(x) mean(x >= 0)),3)

# Fisher skewness of predictive distributions
skdis <- round(sapply(list_pd, function(x) EnvStats::skewness(x)), 3)
skdis[1] <- 0 # symmetric HTS distribution, any asymmetry due to MC noise


## COVID-19 example: Display predictive distributions (Figure 2)
## -----------------------------------------------------------------------------
# Levels and labels of predictive distributions
levels_PD <- c("hts", "nnf", "fixed", "simple", "full")
labels_PD <- c("Higgins-Thompson-Spiegelhalter", "Parametric bootstrap",
               "PCD-fixed", "PCD-simplified", "PCD-full")

# Combine samples from predictive distribution
dat_PD <- data.frame( 
  hts = NA_real_, # display exact density rather than samples
  nnf = pdnnf,
  fixed = pdcovid_fixed$samples[, "theta_new"],
  simple = pdcovid_simplified$samples[, "theta_new"],
  full = pdcovid_full$samples[, "theta_new"]) |>
  pivot_longer(c(hts, nnf, fixed, simple, full), 
               values_to = "s", names_to = "m") |>
  mutate(m = factor(m, levels = c("hts", "nnf", "fixed", "simple", "full")))

# 95% and 99%  equi-tailed and HCD prediction intervals
# * equi-tailed PCD intervals: PredDist internally computes PI's as quantiles; 
#   computed here again for concise code
telePI <- data.frame(
  med = rep(sapply(list_pd, median), times = 2),
  l = c(mapply(function(p) sapply(list_pd, quantile, probs = p), c(0.025, 0.005))),
  u = c(mapply(function(p) sapply(list_pd, quantile, probs = p), c(0.975, 0.995))),
  lev = rep(c("95%", "99%"), each = 5),
  m = factor(rep(levels_PD, times = 2), levels = levels_PD, labels = labels_PD),
  l2 = c(mapply(function(p) sapply(list_pd, function(x) {
    HPDinterval(as.mcmc(x), prob = p)[1]}), c(0.95, 0.99))),
  u2 = c(mapply(function(p) sapply(list_pd, function(x) {
    HPDinterval(as.mcmc(x), prob = p)[2]}), c(0.95, 0.99))))
  
# Determine break points in histogram
brks <- seq(-3.5, 3.5, l = 501)
mids <- (brks[-length(brks)] + brks[-1]) / 2
bin_w <- brks[2] - brks[1]

# Density of HTS predictive distribution
df_hts <- length(es_covid) - 2
mu_hat <- covid_meta$TE.random
se_mu <- sqrt(covid_meta$seTE.random^2 + covid_meta$tau2)
htsdens <- data.frame(
  x = brks,
  f = dt((brks - mu_hat)/se_mu, df = df_hts) / se_mu,
  m = factor("hts", levels = levels_PD, labels = labels_PD),
  region = ifelse(brks < 0, "n", "p")
)

# Plotting data for histograms
plotdat_PD <- dat_PD |>
  mutate(
    m = factor(m, levels = levels_PD, labels = labels_PD),
    bin_idx = as.integer(cut(s, breaks = brks, 
                             include.lowest = TRUE, 
                             right = FALSE))
    ) |>
  group_by(m, bin_idx) |>
  summarise(count = n(), .groups = "drop") |>
  complete(m, bin_idx = seq_len(500), fill = list(count = 0)) |>
  mutate(
    bin_mid   = mids[bin_idx],
    bin_width = bin_w
  ) |>
  group_by(m) |>
  mutate(
    density = count / (sum(count) * bin_width), 
    region  = if_else(bin_mid >= 0, "p", "n")
  ) |> 
  ungroup()

# Blank data.frame to add per panel y-limits
perpanel_ylim <- plotdat_PD |>
  group_by(m) |>
  summarise(x = 0, y = 1.25 * max(density), .groups = "drop")
perpanel_ylim$y[1] <- 2.28 # manual adjustment for HTS

# Y-axis positions of prediction intervals
telePI$ypos <- 0.95 * perpanel_ylim$y

plotdat_PD |>
  ggplot(aes(x = bin_mid, y = density, fill = region, colour = region)) +
  facet_wrap(~ m, nrow = 2, scales = "free") +
  # Histograms
  geom_col(width = bin_w, show.legend = FALSE) +
  # Density and colored area for HTS distribution
  geom_line(data = htsdens, aes(x = x, y = f, color = region)) +
  geom_area(data = htsdens, aes(x = x, y = f, fill = region)) +
  # Color area under predictive distributions
  scale_fill_manual(values = c("n" = "#9FB3C8", "p" = "#4F84B8"), guide = "none") +
  scale_color_manual(values = c("n" = "#9FB3C8", "p" = "#4F84B8"), guide = "none") +
  # New color scale for prediction intervals
  ggnewscale::new_scale_color() +
  # Telescope prediction intervals
  geom_errorbarh(data = telePI |> mutate(type = "HCDP"), inherit.aes = F,
                 aes(y = ypos, xmin = l, xmax = u, alpha = lev, linewidth = lev,
                     color = type), 
                 height = 0) +
  geom_errorbarh(data = telePI |> mutate(type = "Equi-tailed"), inherit.aes = F,
                 aes(y = ypos * 0.925, xmin = l2, xmax = u2, alpha = lev, 
                     linewidth = lev, color = type), 
                 height = 0) +
  # Median of predictive distribution
  geom_point(data = telePI,
             inherit.aes = FALSE, aes(x = med, y = ypos),
             size = 1., color = "#2F5D8A") +
  geom_point(data = telePI,
             inherit.aes = FALSE, aes(x = med, y = ypos * 0.925),
             size = 1., color = "#6FA3C8") +
  # Line to area under density >=0
  geom_segment(
    data = data.frame(
      m = factor(levels_PD, levels = levels_PD, labels = labels_PD),
      x = 1.5, y = 1, 
      xend = c(0.15, 0.2, 0.15, 0.2, 0.2),
      yend = rep(0.1, 5)
    ),
    inherit.aes = FALSE, aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  # Confidence probability of theta_new >= 0
  geom_label(
    data = data.frame(
      m = factor(levels_PD, levels = levels_PD, labels = labels_PD),
      x = 1.5, y = 1, gr0 = gr0),
    inherit.aes = FALSE, aes(x = x, y = y, label = gr0),
    size = 2.5) +
  # Per-panel y-limits with blank data
  geom_blank(data = perpanel_ylim, aes(x = x, y = y), inherit.aes = FALSE) +
  scale_y_continuous(expand = c(0.0, 0.0), breaks = seq(0.0, 3.0, by = 0.5)) +
  scale_alpha_manual(values = c(1, 0.75), guide = "none") +
  scale_x_continuous(limits = c(-2.55, 2.5)) +
  scale_linewidth_manual(values = c(1.0, 0.6)) +
  scale_color_manual(values = c("HCDP" = "#2F5D8A", "Equi-tailed" = "#6FA3C8")) +
  labs(x = expression(theta[new]), y = expression(c(theta[new])),
       color = "Prediction interval", linewidth = "") +
  # Order of legend 
  guides(color = guide_legend(order = 1), linewidth = guide_legend(order = 2)) +
  theme_dk() +
  theme(strip.background = element_rect(fill = "gray99", color = "black"),
        legend.position = c(0.825, 0.19),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(1, "lines"),
        legend.margin = margin(-0.75,0,0,0, unit="cm"),
        legend.text = element_text(size = 9),
        strip.text = element_text(size = 9),
        legend.title = element_text(size = 9, face = "plain"),
        axis.text = element_text(colour = "black", face = "plain",
                                family = "Times", size = 9),
        axis.title = element_text(colour = "black", size = 9, face = "plain"))

## Simulation study results: Load and preprocess results
##------------------------------------------------------------------------------
simresi <- list.files("../results", full.names = TRUE)
dt_v <- do.call(rbind, lapply(simresi, function(i) {
  x <- readRDS(i)
})) |>
  mutate(
    k = factor(k, levels = c(3,5,10,20,50)),
    I2 = factor(I2, levels = c(0, 0.3, 0.6, 0.9),
                labels = c(0, 30, 60, 90)),
  )

# Non-convergences
dt <- dt_v[complete.cases(dt_v), ]
nc_nnf <- dt_v |>
  dplyr::select(starts_with(c("pi.cvr.", "ci.cvr.held.u"))) |>
  summarise(
    hts = mean(is.na(pi.cvr.hts)),
    nnf = mean(is.na(pi.cvr.nnf)),
    fix = mean(is.na(pi.cvr.fix)),
    simple = mean(is.na(pi.cvr.simple)),
    full = mean(is.na(pi.cvr.full)),
    held = mean(is.na(ci.cvr.held.u))
  ) |>
  pull(nnf)

## Sim.-Res.: Coverage of 95% Prediction Intervals (Figure 3)
##------------------------------------------------------------------------------
hlines <- expand.grid(
  k_large = paste0("Large~studies==", 0:2),
  I2 = "iota^2==0*'%'"
)
hlines2 <- expand.grid(
  k_large = paste0("Large~studies==", 0:2),
  I2 = c("iota^2==30*'%'", "iota^2==60*'%'", "iota^2==90*'%'")
)

simplot(dt, "N", "pi.cvr.", nam = 1) +
  geom_hline(
    data = hlines,
    aes(yintercept = 1.00),
    linetype = "dashed", linewidth = 0.25
  ) +
  geom_hline(
    data = hlines2,
    aes(yintercept = 0.95),
    linetype = "dashed", linewidth = 0.25
  ) +
  labs(y = "Coverage of 95% prediction intervals \u00B1 MCSE") +
  scale_y_continuous(labels = scales::label_percent(),
                     breaks = seq(0.8, 1, by = 0.05),
                     expand = c(0.02, 0))

## Sim.-Res.: Pearson correlation between study estimate and PI skewness (Figure 4)
## -----------------------------------------------------------------------------
plotcor(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.hes, nam = 1) +
  labs(y = "Pearson correlation \u00B1 MCSE")

