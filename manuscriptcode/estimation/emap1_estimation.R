## -----------------------------------------------------------------------------
## Code for Manuscript:
##  "Edgington's Method for Random-Effects Meta-Analysis Part I: Estimation"
## by
##   David Kronthaler and Leonhard Held
## -----------------------------------------------------------------------------

## Clear environment
##------------------------------------------------------------------------------
rm(list = ls())

## Library Packages
## -----------------------------------------------------------------------------
library(numDeriv)
library(meta)
library(dplyr)
library(ggplot2) 
library(patchwork)          
library(ggthemes)
library(xtable) 
library(sn)
library(tidyr) 
library(edgemeta)
library(confMeta) 
library(metafor)

## Additional settings
## -----------------------------------------------------------------------------
options(width = 85, digits = 4, show.signif.stars = FALSE)

## Source R Code
##------------------------------------------------------------------------------
# NOTE: set working directory to source file location
source("../functions/utils.R")

## Plotting
##------------------------------------------------------------------------------
# ggplot2 theme
theme_dk <- function() {
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linetype = 1),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(colour = "black", face = "italic",
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
lab_pmu <- expression(p(mu)) # "p(\u03BC)"
lab_mu <- expression(mu) # "\u03BC"
lab_fmu <- expression(c(mu)) # "c(\u03BC)"
lab_tn <- expression(theta[new])
lab_ftn <- expression(f(theta[new]))

# Coloring
clrs <- c("#d31d1d", "#3474b4")

## Example: p-value function and confidence distribution (Glemain et al.,2002)
## -----------------------------------------------------------------------------
ybar <- -0.23
seybar <- 0.59
h0 <- seq(-2.5, 2, l = 250)

# Confidence interval
upex1 <- uniroot(f = function(x) p_fun_wald(x, ybar, seybar) - 0.975,
                          lower = h0[1], upper = h0[250])$root
lpex1 <- uniroot(f = function(x) p_fun_wald(x, ybar, seybar) - 0.025,
                          lower = h0[1], upper = h0[250])$root
# equals 'ybar +- 1.96 * seybar' for the Wald-test

## One-sided p-value function
seg_onesided <- data.frame(
  x = c(ybar, -Inf, upex1, -Inf, lpex1, -Inf),
  xend = c(ybar, ybar, upex1, upex1, lpex1, lpex1),
  y = c(-Inf, 0.5, 0, 0.975, 0, 0.025),
  yend = c(0.5, 0.5, 0.975, 0.975, 0.025, 0.025),
  line_type = factor(c("Median estimate", "Median estimate", 
                       "95% Confidence interval", "95% Confidence interval", 
                       "95% Confidence interval", "95% Confidence interval"),
                     levels = c("Median estimate", "95% Confidence interval"))
)

ponesided <- data.frame(
  mu = h0, 
  p = p_fun_wald(h0, ybar, seybar)) |> 
  ggplot(aes(x = mu, y = p)) +
  geom_vline(xintercept = c(-2, 0, 2), 
             color = "gray90") +
  geom_hline(yintercept = c(0, 0.2, 0.4, 0.6, 0.8, 1),  
             color = "gray90") +
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
    breaks = c(-2, lpex1, ybar, 2, upex1),
    labels = c(
      "-2",
      paste0("<span style='color:#3474b4;'>", 
             round(lpex1, digits = 2), "</span>"),
      paste0("<span style='color:#3474b4;'>", 
             round(ybar, digits = 2), "</span>"),
      "2",
      paste0("<span style='color:#3474b4;'>", 
             round(upex1, digits = 2), "</span>"))) +
  geom_segment(data = seg_onesided, 
               aes(x = x, xend = xend, y = y, yend = yend, linetype = line_type),
               col = clrs[2]) +
  scale_linetype_manual(values = c("Median estimate" = 1, 
                                   "95% Confidence interval" = 2)) +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown(),
        panel.grid.major = element_blank())

## two-sided p-value function
seg_ptwosided <- data.frame(
  x = c(ybar, -Inf, lpex1, upex1, upex1, ybar),
  xend = c(ybar, upex1, lpex1, upex1, upex1, -Inf),
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
  
ptwosided <- data.frame(
  mu = h0, 
  p = p_fun_wald(h0, ybar, seybar, F)
) |>
  ggplot(aes(x = mu, y = p)) +
  geom_vline(xintercept = c(-2, 0, 2), 
             color = "gray90") +
  geom_hline(yintercept = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
             color = "gray90") +
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
  breaks = c(-2, lpex1, ybar, 2, upex1),
  labels = c(
    "-2",
    paste0("<span style='color:#3474b4;'>", 
           round(lpex1, digits = 2), "</span>"),
    paste0("<span style='color:#3474b4;'>", 
           round(ybar, digits = 2), "</span>"),
    "2",
    paste0("<span style='color:#3474b4;'>", 
           round(upex1, digits = 2), "</span>"))) +
   geom_segment(
    data = seg_ptwosided,
    aes(x = x, xend = xend, y = y, yend = yend, linetype = line_type),
    color = clrs[2]
  ) +
  scale_linetype_manual(values = c("Median estimate" = 1, 
                                   "95% Confidence interval" = 2)) +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown(),
        panel.grid.major = element_blank())


## Compute CD values (CDF, density, confidence curve)
pxd <- data.frame(
  h0 = h0,
  cdf = p_fun_wald(h0, ybar, seybar),
  cd = grad(function(x) p_fun_wald(x, ybar, seybar), h0),
  cc = abs(1 - 2 * p_fun_wald(h0, ybar, seybar))) |>
  pivot_longer(c(cdf, cd, cc), values_to = "f", names_to = "d") |>
  mutate(d = factor(d, levels = c("cd", "cdf", "cc")))

# Plot density of CD
pcd2 <- pxd |>
  filter(d == "cd") |>
  ggplot(aes(x = h0, y = f)) + 
  geom_line() +
  labs(x = lab_mu, y = lab_fmu, title = "C") +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dk()

# Plot confidence curve
ci95 <- ybar + c(-1, 1) * qnorm(0.975) * seybar  
ci90 <- ybar + c(-1, 1) * qnorm(0.95)  * seybar  
ci80 <- ybar + c(-1, 1) * qnorm(0.90)  * seybar  
ci_lines <- data.frame(
  x = c(ci95, ci90, ci80),
  yend = rep(c(0.95, 0.9, 0.8), each = 2)
)
ci_lines2 <- data.frame(
  x = rep(-Inf, 3),
  xend = c(ci95[2], ci90[2], ci80[2]),
  y = c(0.95, 0.9, 0.8),
  yend = c(0.95, 0.9, 0.8)
)

pcd3 <- pxd |>
  filter(d == "cc") |>
  ggplot(aes(x = h0, y = f)) + 
  geom_line() +
  geom_segment(data = ci_lines2, inherit.aes = FALSE,
             aes(x = x, xend = xend, y = y, yend = yend),
             linetype = "dashed", color = clrs[2]) +
  geom_segment(data = ci_lines,
               aes(x = x, xend = x, y = 0, yend = yend),
               linetype = "dashed", color = clrs[2]) +
  labs(x = lab_mu, y = expression(CC(mu)), title = "D") +
  scale_y_continuous(expand = c(0.00, 0),
                     breaks = c(0.95, 0.9, 0.8, 0.5, 0.25),
  labels = c("<span style='color:#3474b4;'>0.95</span>",
             "<span style='color:#3474b4;'>0.90</span>",
             "<span style='color:#3474b4;'>0.80</span>",
             "0.5", "0.25")) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_dk() +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown())

# Combine in panel
(ponesided + ptwosided) / (pcd2 + pcd3) + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")


## Example: Serenoa Repens
## -----------------------------------------------------------------------------
load("../data/Serenoa2023.Rdata")

# Fisher's weighted skewness of effect estimates
sereskew <- fwskew(serenoa$effect.size, serenoa$se)

# Random-effects meta-analysis
meme <- metagen(TE = serenoa$effect.size, 
                seTE = serenoa$se, 
                studlab = serenoa$study.name, 
                method.tau = "REML",
                sm = "MD", 
                random = T)

# Edgington's method with additive heterogeneity (by Held et al. 2025)
sere_es <- serenoa$effect.size
sere_se <- sqrt(serenoa$se^2 + meme$tau2)
sere_opt <- remaeffect(es = sere_es, se = sere_se, "NHEU")

# CD-Edgington estimator
me95 <- remaeffect(serenoa$effect.size, 
                   serenoa$se, 
                   level.ci = 0.95, 
                   seed = 982)
conf0 <-  round(mean(me95$cd_mu <= 0), 2) # Confidence probability(mu <= 0)

## Overview Serenoa Data (Table 1)
## -----------------------------------------------------------------------------
serenoa$N <- c(329, 93, 225, 94, 369, 85, 199, 325, 99) # From Franco et al.
sumserenoa <- serenoa |>
  mutate(
    `Estimate` = sprintf("%.2f", effect.size),
    `Standard error` = sprintf("%.2f", se),
    `95% CI` = sprintf("%.2f to %.2f", effect.size - 1.96 * se, effect.size + 1.96 * se)
  ) |>
  select(Study = study.name, N, `Estimate`, `Standard error`, `95% CI`) |>
  xtable(
    caption = "Summary of Serenoa studies analyzed in \\citet{Franco2023}. All studies are 1:1 randomized controlled trials.",
    label = "tab:serenoadata",
    align = c("l", "l", "r", "r", "r", "r"),
    digits = c(0, 0, 0, 2, 2, 0)
  )

# Add footnote
addtorow <- list()
addtorow$pos <- list(nrow(sumserenoa))
addtorow$command <- c(" \n\\multicolumn{5}{l}{\\footnotesize CI = confidence interval.} \\\\ \n")

# Print Table
print(sumserenoa,
      include.rownames = FALSE,
      caption.placement = "top",
      hline.after = NULL,
      booktabs = TRUE,
      add.to.row = addtorow)

## Drapery plot Serenoa
## -----------------------------------------------------------------------------
# ConfMeta meta-analysis with additive heterogeneity
sere_cm <- confMeta(
  serenoa$effect.size,
  serenoa$se,
  heterogeneity = "additive",
  tau2 = meme$tau2,
  study_names = serenoa$study.name,
  conf_level = 0.95,
  fun = p_edgington,
  fun_name = "Edgington  (one-sided input)",
  input_p = "one.sided"
)

# HKSJ method: two-sided p-value function
sere_grid <- seq(-5, 3, l = 1000)
p2hk <- unlist(lapply(sere_grid, function(x) {
  s2 <- metagen(serenoa$effect.size, 
                serenoa$se, 
                method.tau = "REML", 
               method.random.ci = "HK", 
               null.effect = x)$pval.random
  s1 <- if (s2 <= 0.5) s2 / 2 else 1 - s2 / 2
  return(2 * min(s1, 1 - s1))
}))

# CD-Edgington: two-sided p-value function
CDFscd <- ecdf(me95$cd_mu)
p2scd <- ifelse(CDFscd(sere_grid) <= 0.5,
                CDFscd(sere_grid) * 2, 
                2 * (1 - CDFscd(sere_grid)))

# 95% and 99% confidence intervals for all methods
LCI <- c(0.99, 0.95)

# Edgington with additive heterogeneity
CIedgington <- lapply(LCI, function(x) { 
  remaeffect(sere_es, serenoa$se, "NHEU", level.ci = x) 
})

# Classical REMA
CIre <- lapply(LCI * 100, function(x) {
  rma(yi = sere_es, sei = serenoa$se, method = "REML", level = x)
})

# Monte Carlo CD-Edgington
CIsampled <- lapply(LCI, function(lvl) {
  remaeffect(sere_es, serenoa$se, level.ci = lvl, seed = 982)
})

# GAQ CD-Edgington
CIaq <- remaeffect(sere_es, serenoa$se, method = "GAQ")

# Hartung-Knapp
CIhk <- lapply(LCI * 100, function(lvl) {
  rma(yi = sere_es, sei = serenoa$se, method = "REML", level = lvl, test = "knha")
})

# Telescope lines data
CIdt <- data.frame(
  Method = factor(rep(c("Random-effects", 
                        "Hartung-Knapp-Sidik-Jonkman",
                        "Edgington", 
                        "CD-Edgington"),
                      each = 2),
                  levels = c("CD-Edgington", 
                             "Edgington", 
                             "Hartung-Knapp-Sidik-Jonkman", 
                             "Random-effects")),
  cl = factor(rep(c("99%", "95%"), times = 4), levels = c("99%", "95%")),
  p = c(
    rep(meme$TE.random, 4),                   
    rep(sere_opt$estimate, 2),                   
    rep(CIsampled[[2]]$estimate, 2)            
  ),
  l = c(
    sapply(CIre, function(x) x$ci.lb),
    sapply(CIhk, function(x) x$ci.lb),
    sapply(CIedgington, function(x) x$CI[1]),
    sapply(CIsampled, function(x) x$CI[1])
  ),
  u = c(
    sapply(CIre, function(x) x$ci.ub),
    sapply(CIhk, function(x) x$ci.ub),
    sapply(CIedgington, function(x) x$CI[2]),
    sapply(CIsampled, function(x) x$CI[2])
  )
)
CIdt$ypos <- rep(seq(1.02,1.12, l = 4), each = 2)

# Coloring
cdrap <- c("Random-effects" = "#6C91BF", 
           "Hartung-Knapp-Sidik-Jonkman" = "#90A4AE", 
           "Edgington" = "#88C0A9", 
           "CD-Edgington" = "#F1C06E")

# Drapery plot
drape <- autoplot(sere_cm)[[1]]
drapeprint <- layer_data(drape, 2) |>
  ggplot(aes(x = x, y = y, group = group)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "solid") +
  geom_hline(yintercept = 0.05, color = "black", linewidth = 0.75, linetype = "dashed") +
  # Study p-value functions
  geom_line(linetype = "dashed", color = "lightgrey", linewidth = 0.5) +
  # Classical REMA
  geom_line(data = layer_data(drape, 3), aes(x = x, y = y, group = group), 
            color = cdrap["Random-effects"], linewidth = 0.75, linetype = 1) +
  # Edgington
  geom_line(data = layer_data(drape, 7), aes(x = x, y = y, group = group), 
            color = cdrap["Edgington"], linewidth = 0.75, linetype = 1) +
  # Hartung-Knapp
  geom_line(inherit.aes = FALSE, linewidth = 1, linetype = 1,
            data = data.frame(x = sere_grid, y = p2hk, 
                              Method = "Hartung-Knapp-Sidik-Jonkman"),
            aes(x = x, y = y, color = Method)) +
  # CD-Edgington
  geom_line(inherit.aes = FALSE, linewidth = 1,
            data = data.frame(x = sere_grid, y = p2scd, 
                              Method = "CD-Edgington"),
            aes(x = x, y = y, color = Method)) +
  labs(x = lab_mu) +
  coord_cartesian(xlim = c(-3, 1), ylim = c(0, 1.15)) +
  # Dummy layers for legend
  geom_line(data = data.frame(x = -10, y = -10, 
                              Method = "Random-effects"),
            aes(x = x, y = y, color = Method, group = Method)) +
  geom_line(data = data.frame(x = -10, y = -10,
                              Method = "Edgington"),
            aes(x = x, y = y, color = Method, group = Method)) +
  scale_color_manual(name = "Method", values = cdrap) +
  # Confidence intervals
  geom_errorbarh(data = CIdt, aes(y = ypos, x = p, xmin = l, xmax = u, 
                                  group = Method, linewidth = cl, alpha = cl,
                                  color = Method), height = 0) +
  geom_point(data = CIdt, aes(y = ypos, x = p, group = Method, color = Method),
             size = 2) +
  scale_linewidth_manual(values = c(0.75, 1.25)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_y_continuous(
    name = lab_pmu,
    sec.axis = sec_axis( trans=~.*-100 + 100, name="Confidence level [%]",
                         breaks = c(0, 20, 40, 60, 80, 95)),
    breaks = c(0.05, 0.2, 0.4, 0.6, 0.8, 1.0),
    expand = c(0,0)) +
  theme_dk() +
  theme(axis.title.y.right = element_text(angle = 90))

drapeprint

## Estimated confidence distributions (tau2 and mu)
## -----------------------------------------------------------------------------
ctau2 <- data.frame(x = CIsampled[[2]]$cd_tau2,
           sere_grid = seq(0, 10, l = 1000),
           f = cdtau2( seq(0, 10, l = 1000), serenoa$effect.size, serenoa$se)) |>
  ggplot(aes(x = x)) + 
  geom_histogram(bins = 350, aes(y = ..density..), fill = "gray95", color = "gray20") +
  geom_line(inherit.aes = FALSE, mapping = aes(x = sere_grid, y = f), color = "royalblue",
            linewidth = 1) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(y = expression(c~(tau^2)), x = expression(tau^2), title = "A") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.01)) +
  theme_dk()

hdt <- hist(CIsampled[[2]]$cd_mu, breaks = 150, plot = FALSE)
dfh <- data.frame(
  mid = hdt$mids,
  density = hdt$density,
  filler = ifelse(hdt$mids < 0, "lightblue", "gray95")
)

dfl <- data.frame(
  sere_grid = seq(-3, 1.5, l = 1000),
  f = CIaq$fcd(seq(-3, 1.5, l = 1000))
)

cmu <- ggplot(dfh, aes(x = mid, y = density)) +
  geom_col(aes(fill = filler), color = "gray20", width = diff(hdt$breaks)[1]) +
  geom_line(data = dfl, aes(x = sere_grid, y = f), color = "royalblue", linewidth = 1) +
  coord_cartesian(xlim = c(-3,1.5)) +
  labs(y = lab_fmu, x = lab_mu, title = "B") +
  geom_segment(x = -2.1, y = 0.75, xend = -1, yend = 0.25,
               color = "black") +
  geom_label(
    aes(label = paste0("Conf(mu < 0) == ", conf0)),
    x = -2.1, y = 0.75,
    parse = TRUE,
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.01)) +
  theme_dk()

ctau2 + cmu

## Point estimates and confidence intervals for average effect (Table 2)
## -----------------------------------------------------------------------------
# Hartung-Knapp method
sere_meta_hk <- metagen(sere_es, serenoa$se, random = TRUE, method.random.ci = "HK")

# Bind together and display
estiserenoa <- CIdt |>
  filter(cl == "95%") |>
  bind_rows(data.frame(
    Method = "CD-Edgington (GAQ)",
    cl = "95%",
    p = CIaq$estimate,
    l = CIaq$CI[1],
    u = CIaq$CI[2]
  )) |>
  rowwise() |>
  mutate(
    Method = ifelse(Method == "CD-Edgington", "CD-Edgington (MC)", Method),
    sk = skew_pi(c(l, u), p),
    ci = paste0(sprintf("%.2f",l), "  to  ", sprintf("%.2f",u)),
    w = u - l
  ) |>
  bind_cols(pval = biostatUZH::formatPval(c(
      metagen(serenoa$effect.size, serenoa$se)$pval.random,
      metagen(serenoa$effect.size, serenoa$se, method.random.ci = "HK")$pval.random,
      sere_cm$p_0[2],
      CIsampled[[2]]$pval,
      CIaq$pval
    ))) |>
  dplyr::select(Method, p, ci, w, sk, pval) |>
  rename(
    "Estimator" = "Method",
    "Estimate" = "p", 
    "95\\% CI" = "ci",
    "Skewness" = "sk",
    "Width" = "w",
    "p-value" = "pval"
  ) |> 
  xtable(caption = "Point estimates and 95\\% confidence intervals for the average treatment effect across Serenoa studies. $P$-values are computed for the null hypothesis $\\mu_0 = 0$.",
         align = c("l", "l", "r", "r", "r", "r", "r"),
         label = "tab:serepointestimates")

# Add footnote
addtorow <- list()
addtorow$pos <- list(nrow(estiserenoa))
addtorow$command <- c(" \n\\multicolumn{6}{l}{\\footnotesize CI = confidence interval, GAQ = global adaptive quadrature, MC = Monte Carlo.} \\\\ \n")

# Print table
print(estiserenoa,
      include.rownames = FALSE,
      caption.placement = "top",
      hline.after = NULL,
      booktabs = TRUE,
      sanitize.text.function = identity,
      add.to.row = addtorow)

## Simulation results: Load and preprocess results
##------------------------------------------------------------------------------
simresi <- list.files("../../results")
dt <- do.call(rbind, lapply(simresi, function(i) {
  x <- readRDS(paste0("../../results/", i))
})) |>
  mutate(
    k = factor(k, levels = c(3,5,10,20,50)),
    I2 = factor(I2, levels = c(0, 0.3, 0.6, 0.9),
                labels = c(0, 30, 60, 90)), 
      across(
      starts_with("q05"),
      ~ ifelse(I2 == 0, (. - 0)^2, (. - 0.5)^2),
      .names = "{.col}"
    ),
    across(
      starts_with("q01"),
      ~ ifelse(I2 == 0, (. - 0)^2, (. - 0.1)^2),
      .names = "{.col}"
    ),
    across(
      starts_with("q005"),
      ~ ifelse(I2 == 0, (. - 0)^2, (. - 0.05)^2),
      .names = "{.col}"
    )
  ) |>
  rename("sqe.held.a" = "seq.held.a")

## ReSim: Coverage 95% confidence intervals
##------------------------------------------------------------------------------
simplot(dt, "N", "ci.cvr.", nam = 4) +
  labs(y = "Coverage of 95% confidence intervals \u00B1 MCSE") + 
  geom_hline(yintercept = 0.95, linetype = "dashed", size = 0.25) +
  scale_y_continuous(labels = scales::label_percent(),
                     breaks = c(0.7, 0.8, 0.9, 1))

## ReSim: Width of 95% confidence intervals
## -----------------------------------------------------------------------------
simplot(dt = dt, distr = "N", stw_var = "ci.w.", nam = 4) +
  labs(y = "Width of 95% confidence intervals \u00B1 MCSE")

## Maximum avg. bias of point estimators
## -----------------------------------------------------------------------------
# Normal effects
bmax <- min(dt |>
              filter(dist == "N") |>
              group_by(I2, k, k_large) |>
              summarise(
                bmax = mean(bias.ivw.random)
              ) |>
              pull(bmax))
# Skew normal effects
bmaxlsn <- max(dt |>
                 filter(dist == "LSN") |>
                 group_by(I2, k, k_large) |>
                 summarise(
                   bmax = mean(bias.held.a)
                 ) |>
                 pull(bmax))

## ReSim: Bias point estimator mu
##------------------------------------------------------------------------------
simplot(dt = dt, distr = "N", stw_var = "bias.", nam = 3) +
  labs(y = "Bias \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  scale_y_continuous(limits = c(-0.05, 0.05))


################################################################################
################################################################################


## -----------------------------------------------------------------------------
## Code for Supplementary Material of Manuscript:
##  "Edgington's Method for Random-Effects Meta-Analysis Part I: Estimation"
## by
##   David Kronthaler and Leonhard Held
## -----------------------------------------------------------------------------

## Serenoa Repens and Corticosteroids & Covid19 Examples: CD of tau2
##------------------------------------------------------------------------------
# Serenoa Repens
sere_mameta <- metagen(serenoa$effect.size, 
                       serenoa$se,
                       random = TRUE,
                       method.tau2 = "REML")

# Corticosteroids and Covid19
load("../data/HeldPawelHofmann2025.Rdata")
hph_ftau2_tau2 <- metagen(datHPH2025$logOR, 
                          datHPH2025$logSE, 
                          random = TRUE,
                          method.tau2 = "REML")

## CD Serenoa Repens
## -----------------------------------------------------------------------------
sxi_ftau2_sere <- seq(0, 5, l = 1000)

# Compute (and scale to account for integration error) CDF
Ftau2_sere <- sapply(sxi_ftau2_sere, function(xi) {
  result <- integrate( function(x) {
    cdtau2(x, serenoa$effect.size, serenoa$se)},
    lower = 0, upper = xi,
    rel.tol = 1e-8, abs.tol = 1e-12)
  result$value
})
Ftau2_sere <- Ftau2_sere / max(Ftau2_sere)

# Quantile function
QuantTau2_sere <- approxfun(Ftau2_sere, sxi_ftau2_sere, rule = 2)

# Compute median and 95% confidence interval
median_q <- QuantTau2_sere(0.5)
q025 <- QuantTau2_sere(0.025)
q975 <- QuantTau2_sere(0.975)

# Compute and plot confidence density
ftau2_sere <- data.frame(
  xi = sxi_ftau2_sere,
  x2 = cdtau2(sxi_ftau2_sere, serenoa$effect.size, serenoa$se)) |>
  ggplot(aes(x = xi, y = x2)) +
  geom_segment(aes(x = median_q, xend = median_q, y = -Inf, yend = Inf),
               size = 0.3, linetype = 1, color = clrs[2]) +
  geom_segment(aes(x = q025, xend = q025, y = -Inf, yend = Inf),
               size = 0.3, linetype = 2, color = clrs[2]) +
  geom_segment(aes(x = q975, xend = q975, y = -Inf, yend = Inf),
               size = 0.3, linetype = 2, color = clrs[2]) +
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
                      format(q975, digits = 2), "</span>"),
               "4", "5")) +
  labs(
    y = expression(c~(tau^2)), x = expression(tau^2),
    title = "A"
  ) +
  theme_dk() +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown())

# plot CDF
pFtau2_Sere <- data.frame(x = sxi_ftau2_sere, cdf = Ftau2_sere) |>
  ggplot(aes(x = x, y = cdf)) +
  geom_line() +
  geom_vline(xintercept = c(0, 1, 2, 3, 4, 5),
             color = "gray90") +
  geom_hline(yintercept = c(0, 0.25, 0.5, 0.75, 1),
             color = "gray90") +
  scale_y_continuous(expand = c(0.005,0),
                     breaks = c(0.025, 0.25, 0.5, 0.75, 0.975),
                     labels = c("<span style='color:#3474b4;'>0.025</span>", "0.25",
                                "<span style='color:#3474b4;'>0.500</span>", "0.75",
                                "<span style='color:#3474b4;'>0.975</span>")) +
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
                      format(q975, digits = 2), "</span>"),
               "4", "5")) +
  geom_segment(aes(y = 0.5, yend = 0.5, x = 0, xend = median_q,
                   linetype = "Median estimate"),
               color = clrs[2], size = 0.3) +
  geom_segment(aes(x = median_q, xend = median_q, y = 0, yend = 0.5,
                   linetype = "Median estimate"),
               color = clrs[2], size = 0.3) +
  geom_segment(aes(y = 0.025, yend = 0.025, x = 0, xend = q025,
                   linetype = "95% Confidence interval"),
               color = clrs[2], size = 0.3) +
  geom_segment(aes(y = 0.975, yend = 0.975, x = 0, xend = q975,
                   linetype = "95% Confidence interval"),
               color = clrs[2], size = 0.3) +
  geom_segment(aes(x = q025, xend = q025, y = 0, yend = 0.025,
                   linetype = "95% Confidence interval"),
               color = clrs[2], size = 0.3) +
  geom_segment(aes(x = q975, xend = q975, y = 0, yend = 0.975,
                   linetype = "95% Confidence interval"),
               color = clrs[2], size = 0.3) +
  scale_linetype_manual(values = c(
    "Median estimate" = 1,
    "95% Confidence interval" = 2
  ), guide = guide_legend(reverse = TRUE)) +
  labs(
    y = expression(C~(tau^2)), x = expression(tau^2),
    title = "C"
  ) +
  theme_dk() +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown(),
        panel.grid.major = element_blank())

## CD Corticosteroids and Covid19
## -----------------------------------------------------------------------------
sxi_ftau2_hph <- seq(0, 2, l = 1000)

# Compute (and scale to account for integration error) CDF
Ftau2_hph <- sapply(sxi_ftau2_hph, function(xi) {
  result <- integrate( function(x) {
    cdtau2(x, datHPH2025$logOR, datHPH2025$logSE)},
    lower = 0, upper = xi,
    rel.tol = 1e-8, abs.tol = 1e-12)
  result$value
})
Ftau2_hph <- Ftau2_hph / max(Ftau2_hph)

# Compute median and 95% confidence interval
QuantTau2_hph <- approxfun(Ftau2_hph, sxi_ftau2_hph, rule = 2)
median_q_hph <- QuantTau2_hph(0.5)
q025_hph <- QuantTau2_hph(0.025)
q975_hph <- QuantTau2_hph(0.975)

# plot confidence density
ftau2_hph <- data.frame(
  xi = sxi_ftau2_hph,
  x2 = cdtau2(sxi_ftau2_hph, datHPH2025$logOR, datHPH2025$logSE)) |>
  ggplot(aes(x = xi, y = x2)) +
  geom_line() +
  geom_segment(aes(x = median_q_hph, xend = median_q_hph, y = -Inf, yend = Inf),
               size = 0.3, linetype = 1, color = clrs[2]) +
  geom_segment(aes(x = q025_hph, xend = q025_hph, y = -Inf, yend = Inf),
               size = 0.3, linetype = 2, color = clrs[2]) +
  geom_segment(aes(x = q975_hph, xend = q975_hph, y = -Inf, yend = Inf),
               size = 0.3, linetype = 2, color = clrs[2]) +
  scale_y_continuous(expand = c(0.005,0)) +
  scale_x_continuous(
    expand = c(0.005, 0),
    breaks = c(q025_hph, median_q_hph, 0.5, 1, 1.5, q975_hph, 2),
    labels = c(
      paste0("<span style='color:#3474b4;'>", format(q025_hph, digits = 2),
             "</span>"),
      paste0("<span style='color:#3474b4;'><br>",
             format(median_q_hph, digits = 2),
             "</span>"),
      "0.5", "1", "1.5",
      paste0("<span style='color:#3474b4;'><br>",
             format(q975_hph, digits = 2),
             "</span>"),
      "2")) +
  labs(
    y = expression(c~(tau^2)), x = expression(tau^2),
    title = "B"
  ) +
  theme_dk() +
  theme(legend.position = "none")  +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown())

# plot CDF
pFtau2_hph <- data.frame(x = sxi_ftau2_hph, cdf = Ftau2_hph) |>
  ggplot(aes(x = x, y = cdf)) +
  geom_vline(xintercept = c(0, 0.5, 1, 1.5, 2),
             color = "gray90") +
  geom_hline(yintercept = c(0, 0.25, 0.5, 0.75, 1),
             color = "gray90") +
  geom_line() +
  scale_y_continuous(expand = c(0.005,0),
                     breaks = c(0.025, 0.25, 0.5, 0.75, 0.975),
                     labels = c("<span style='color:#3474b4;'>0.025</span>", "0.25",
                                "<span style='color:#3474b4;'>0.500</span>", "0.75",
                                "<span style='color:#3474b4;'>0.975</span>")) +
  scale_x_continuous(
    expand = c(0.005, 0),
    breaks = c(q025_hph, median_q_hph, 0.5, 1, 1.5, q975_hph, 2),
    labels = c(
      paste0("<span style='color:#3474b4;'>", format(q025_hph, digits = 2),
             "</span>"),
      paste0("<span style='color:#3474b4;'><br>",
             format(median_q_hph, digits = 2),
             "</span>"),
      "0.5", "1", "1.5",
      paste0("<span style='color:#3474b4;'><br>",
             format(q975_hph, digits = 2),
             "</span>"),
      "2")) +
  geom_segment(aes(y = 0.5, yend = 0.5, x = 0, xend = median_q_hph,
                   linetype = "Median estimate"), color = clrs[2], size = 0.3) +
  geom_segment(aes(x = median_q_hph, xend = median_q_hph, y = 0, yend = 0.5,
                   linetype = "Median estimate"), color = clrs[2], size = 0.3) +
  geom_segment(aes(y = 0.025, yend = 0.025, x = 0, xend = q025_hph,
                   linetype = "95% Confidence interval"),
               color = clrs[2], size = 0.3) +
  geom_segment(aes(y = 0.975, yend = 0.975, x = 0, xend = q975_hph,
                   linetype = "95% Confidence interval"),
               color = clrs[2], size = 0.3) +
  geom_segment(aes(x = q025_hph, xend = q025_hph, y = 0, yend = 0.025,
                   linetype = "95% Confidence interval"),
               color = clrs[2], size = 0.3) +
  geom_segment(aes(x = q975_hph, xend = q975_hph, y = 0, yend = 0.975,
                   linetype = "95% Confidence interval"),
               color = clrs[2], size = 0.3) +
  scale_linetype_manual(values = c(
    "Median estimate" = 1,
    "95% Confidence interval" = 2
  ), guide = guide_legend(reverse = TRUE)) +
  labs(
    y = expression(C~(tau^2)), x = expression(tau^2),
    title = "D"
  ) +
  theme_dk() +
  theme(axis.text.y = ggtext::element_markdown(),
        axis.text.x = ggtext::element_markdown(),
        panel.grid.major = element_blank())

# Align plots
ftau2_sere + ftau2_hph + pFtau2_Sere + pFtau2_hph +
  plot_layout(axis_titles = "collect_x", guides = "collect") &
  theme(legend.position = "bottom",
        legend.text = element_text(size = 10)) &
  guides()
# Warnings can be ignored here

## Pilot simulation comparing Monte Carlo and GAQ CD-Edgington estimators
## -----------------------------------------------------------------------------
if (!file.exists("../data/MCvsGAQ.RDS")) {
  library(doParallel)
  library(foreach)
  
  # One simulation iteration
  os <- function(k, I2) {
    
    # Generate data
    ni <- 50
    mu <- -0.3
    se <- sqrt(rchisq(k, df = 2 * (ni - 1)) * ((ni - 1) * ni)^(-1))
    tau2 <- 1 / k * sum(2 / ni) * (I2 / (1 - I2))
    es <- rnorm(k, mean = mu, sd = sqrt(tau2))
    hes <- rnorm(k, es, sqrt(2/ni)) 
    
    # Monte Carlo
    mc <- remaeffect(hes, se, method = "MC")
    
    # GAQ
    aq <-  remaeffect(hes, se, method = "GAQ")
    
    # Return
    return(unname(c(mc$estimate, mc$CI, aq$estimate, aq$CI)))
  }
  
  # Number of iterations
  niter <- 1000
  
  # Grid of investigated conditions
  presimlev <- expand.grid(
    k = c(3, 5, 10, 20, 50),
    I2 = c(0.0, 0.3, 0.6, 0.9),
    i = 1:niter
  )
  
  # Parallelization
  num_cores <- parallel::detectCores() - 5
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  respresim <- foreach(j = seq_len(nrow(presimlev)), 
                       .errorhandling = "remove",
                       .packages = "metaprediction",
                       .combine = rbind) %dopar% {
                         os(presimlev$k[j], presimlev$I2[j])
                       }
  stopCluster(cl)
  
  # Result
  respresim <- cbind(respresim, presimlev)
  colnames(respresim) <- c("mcp", "mcl", "mcu", "aqp", "aql", "aqu",
                           "k", "I2", "i")
  saveRDS(respresim, file = "../data/MCvsGAQ.RDS")
} else {
  respresim <- readRDS("../data/MCvsGAQ.RDS")  
}   

## Supplementary Table 1: Mean diffs. in estimates and CI limits between MC and GAQ
## -----------------------------------------------------------------------------
tt <- respresim %>%
  mutate(
    mc_cover = mcl <= -0.3 & mcu >= -0.3,
    aq_cover = aql <= -0.3 & aqu >= -0.3,
    mc_bias = mcp + 0.3,
    aq_bias = aqp + 0.3
  ) |> 
  group_by(k, I2) %>% 
  summarise(
    mcbias  = mean(mc_bias),  mc_bias_se  = sd(mc_bias) / sqrt(n()),
    aqbias = mean(aq_bias), aq_bias_se = sd(aq_bias) / sqrt(n()),
    mccover = mean(mc_cover), mc_cover_se = sd(mc_cover) / sqrt(n()),
    aqcover = mean(aq_cover), aq_cover_se = sd(aq_cover) / sqrt(n()),
    diffp = mean(mcp - aqp), diff_p_se = sd(mcp - aqp) / sqrt(n()),
    difflower = mean(mcl - aql), diff_lower_se = sd(mcl - aql)/sqrt(n()),
    diffupper = mean(mcu - aqu), diff_upper_se = sd(mcu - aqu)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    mcbias  = sprintf("%.3f [%.3f]", mcbias,  mc_bias_se),
    aqbias = sprintf("%.3f [%.3f]", aqbias, aq_bias_se),
    mccover = sprintf("%.3f [%.3f]", mccover, mc_cover_se),
    aqcover = sprintf("%.3f [%.3f]", aqcover, aq_cover_se),
    diffp = sprintf("%.3f [%.3f]", diffp, diff_p_se),
    difflower = sprintf("%.3f [%.3f]", difflower, diff_lower_se),
    diffupper = sprintf("%.3f [%.3f]", diffupper, diff_upper_se),
  ) %>%
  select(k, I2, mcbias, aqbias, mccover, aqcover, diffp, difflower, diffupper)

fm <- function(dt, metric) {
  dt %>%
    select(k, I2, !!sym(metric)) %>%
    pivot_wider(names_from = I2, values_from = !!sym(metric)) %>%
    arrange(k) %>%
    mutate(k = as.character(k)) %>%
    rename_with(~ paste0(as.numeric(.x) * 100, "%"), -k)
}

header <- function(dt, section_name) {
  hr <- as.list(c(section_name, "", "", "", ""))
  i2r <- as.list(c("Studies ($k$)", "$I^2$ = 0\\%", "30\\%", "60\\%", "90\\%"))
  rbind(hr, i2r, dt)
}

MCGAQ <- xtable(bind_rows(
  header(fm(tt, "diffp"),  "Estimate (MC - GAQ)"),
  header(fm(tt, "difflower"), "Lower 95\\% CI (MC - GAQ)"),
  header(fm(tt, "diffupper"), "Upper 95\\% CI (MC - GAQ)")),
  caption = "Mean differences with Monte Carlo standard errors in point estimates and 95\\% confidence interval limits between Monte Carlo sampling and global adaptive quadrature integretation approaches.", 
  label = "tab:MCvsGAQ")

# Add footnote
addtorow <- list()
addtorow$pos <- list(nrow(MCGAQ))
addtorow$command <- c(" \n\\multicolumn{5}{l}{\\footnotesize CI = confidence interval, GAQ = global adaptive quadrature, MC = Monte Carlo.} \\\\ \n")

# Print table
print(MCGAQ,
      include.rownames = FALSE, 
      include.colnames = FALSE,
      caption.placement = "top",
      hline.after = NULL, 
      booktabs = TRUE,
      sanitize.text.function = identity,
      add.to.row = addtorow)

## Supplementary Table 2: Bias and coverage of MC and GAQ
## -----------------------------------------------------------------------------
MCGAQ_bc <- xtable(bind_rows(
  header(fm(tt, "mcbias"),  "MC: Bias"), 
  header(fm(tt, "aqbias"), "GAQ: Bias"),
  header(fm(tt, "mccover"), "MC: 95\\% CI coverage"),
  header(fm(tt, "aqcover"), "GAQ: 95\\% CI coverage")),
  caption = "Bias of point estimators and coverage of 95\\% confidenc intervals with Monte Carlo standard errors for Monte Carlo sampling and global adaptive quadrature integretation approaches.", 
  label = "tab:MCvsGAQ_bc")  

# Add footnote
addtorow <- list()
addtorow$pos <- list(nrow(MCGAQ_bc))
addtorow$command <- c(" \n\\multicolumn{5}{l}{\\footnotesize CI = confidence interval, GAQ = global adaptive quadrature, MC = Monte Carlo.} \\\\ \n")

# Print table
print(MCGAQ_bc, 
      include.rownames = FALSE,
      include.colnames = FALSE,
      caption.placement = "top",
      hline.after = NULL, 
      booktabs = TRUE,
      sanitize.text.function = identity,
      add.to.row = addtorow)


## Example of skew-normal distribution
##------------------------------------------------------------------------------
# Parameters by moment-matching
mu = -0.3
tau2 = 0.5
alpha <- - 4
d <- alpha / sqrt(1 + alpha ^ 2)
omega <- sqrt(tau2 / (1 - 2 * (d^2) / pi))
xi <- mu - omega * d * sqrt(2 / pi)

# Parameters, grid and evaluate density
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

# Plot
ggplot(pars_sn, aes(x = x, y = dens, color = alpha, group = alpha)) +
  facet_grid(~label,
             labeller = label_parsed) +
  geom_line(color = "black", size = 0.5) +
  labs(y = "Density", x = "Y") +
  scale_y_continuous(expand = c(0,0.01), limits = c(0,0.65)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_dk() +
  theme(
    strip.background = element_rect(fill = "gray99", color = "black"),
    legend.position = "none"
  )

## ReSim: Coverage 95% confidence intervals
##------------------------------------------------------------------------------
# Skew-normal effects
simplot(dt, "LSN", "ci.cvr.", nam = 4) +
  labs(y = "Coverage of 95% confidence intervals \u00B1 MCSE") +
  geom_hline(yintercept = 0.95, linetype = "dashed", size = 0.25) +
  scale_y_continuous(labels = scales::label_percent())

## ReSim: Width of 95% confidence intervals
## -----------------------------------------------------------------------------
# Skew-normal effects
simplot(dt = dt, distr = "LSN", stw_var = "ci.w.", nam = 4) +
  labs(y = "Width of 95% confidence intervals \u00B1 MCSE")

## ReSim: Skewness agreement of CIs with effect estimates
## -----------------------------------------------------------------------------
# Normal effects
plotcor(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.hes, nam = 2) +
  labs(y = "Pearson correlation \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning here is due to the classical approach always being symmetric, and 
# correlations of skewness yielding NA.

## ReSim: Skewness agreement of CIs with effect estimates
## -----------------------------------------------------------------------------
# Skew-normal effects
plotcor(dt = dt, distr = "LSN", stw_var = "ci.sk.", var = sk.hes, nam = 2) +
  labs(y = "Pearson correlation \u00B1 MCSE")  +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning here is due to the classical approach always being symmetric, and 
# correlations of skewness yielding NA.

## ReSim: Skewness sign agreement of CIs with effect estimates (Cohens kappa)
## -----------------------------------------------------------------------------
# Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.hes, nam = 2) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning here is due to the classical approach always being symmetric, and 
# correlations of skewness yielding NA.

## ReSim: Skewness sign agreement of CIs with effect estimates (Cohens Kappa)
## -----------------------------------------------------------------------------
# Skew-Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.hes, nam = 2) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning here is due to the classical approach always being symmetric, and 
# correlations of skewness yielding NA.

## ReSim: Skewness agreement confidence intervals with true effects
## -----------------------------------------------------------------------------
# Normal effects
plotcor(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.es, nam = 2) +
  labs(y = "Pearson correlation \u00B1 MCSE")  +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) 
# Warning here is due to the classical approach always being symmetric, and 
# correlations of skewness yielding NA.

## ReSim: Skewness agreement confidence intervals with true effects
## -----------------------------------------------------------------------------
# Skew-normal effects
plotcor(dt = dt, distr = "LSN", stw_var = "ci.sk.", var = sk.es, nam = 2) +
  labs(y = "Pearson correlation \u00B1 MCSE")  +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) 
# Warning here is due to the classical approach always being symmetric, and 
# correlations of skewness yielding NA.

## ReSim: Skewness sign agreement of CIs with true effects (Cohens Kappa)
## -----------------------------------------------------------------------------
# Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.es, nam = 2) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning here is due to the classical approach always being symmetric, and 
# correlations of skewness yielding NA.

## ReSim: Skewness sign agreement of CIs with true effects (Cohens Kappa)
## -----------------------------------------------------------------------------
# Skew-Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "ci.sk.", var = sk.es, nam = 2) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning here is due to the classical approach always being symmetric, and 
# correlations of skewness yielding NA.

## ReSim: Bias point estimator mu
##------------------------------------------------------------------------------
# Skew-normal effects
simplot(dt = dt, distr = "LSN", stw_var = "bias.", nam = 3) +
  labs(y = "Bias \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## ReSim: MSE point estimator mu
##------------------------------------------------------------------------------
# Normal effects
simplot(dt, "N", "sqe.", nam = 3) +
  labs(y = "MSE \u00B1 MCSE")+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## ReSim: MSE point estimator mu
##------------------------------------------------------------------------------
# Skew-normal effects
simplot(dt, "LSN", "sqe.", nam = 3) +
  labs(y = "MSE \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

