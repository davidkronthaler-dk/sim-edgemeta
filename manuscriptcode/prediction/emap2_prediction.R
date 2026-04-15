## -----------------------------------------------------------------------------
## Code for Manuscript:
##  "Edgington's Method for Random-Effects Meta-Analysis Part II: Prediction"
## by
##   David Kronthaler and Leonhard Held
## -----------------------------------------------------------------------------

## Clear environment
##------------------------------------------------------------------------------
rm(list = ls())
 
## Library Packages
## -----------------------------------------------------------------------------
library(meta)
library(dplyr)
library(ggplot2)
library(patchwork)          
library(ggthemes) 
library(xtable) 
library(sn)
library(latex2exp)
library(tidyr)  
library(Rcpp)
library(edgemeta)
library(confMeta) 
library(metafor)
library(coda)

## Additional settings 
## -----------------------------------------------------------------------------
options(width = 85, digits = 4, show.signif.stars = FALSE)

## Source R Code
##------------------------------------------------------------------------------
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
lab_pmu <- "p(\u03BC)"
lab_mu <- "\u03BC"
lab_fmu <- "c(\u03BC)"
lab_tn <- expression(theta[new])
lab_ftn <- expression(c(theta[new]))

## COVID19 Example: Load Data, Classical REMA
##------------------------------------------------------------------------------
# Load data
load("../data/HeldPawelHofmann2025.Rdata")
escov <- datHPH2025$logOR  # extract study effects
secov <- datHPH2025$logSE  # extract study standard errors

# Fisher's weighted skewness of effect estimates
covskew <- round(fwskew(escov, secov), 3)

# REMA with HTS prediction interval
ma_meta_sere <- metagen(TE = escov,
                        seTE = secov,
                        method.tau = "PM", 
                        method.random.ci = "HK",
                        common = F,
                        method.predict = "HTS")

# CD-Edgington estimator
me95 <- remaeffect(escov, secov, level.ci = 0.95, seed = 982)

# Skewness of Edgington's confidence distribution
skECD <- EnvStats::skewness(me95$cd_mu)

# Confidence of avg. effect < 0
confcov <- mean(me95$cd_mu < 0)

# Estimates of average effect for text
EShk <- paste0(round(ma_meta_sere$TE.random, 2),
               " (95\\% confidence interval from ", 
               round(ma_meta_sere$lower.random, 2), " to ", 
               round(ma_meta_sere$upper.random, 3),
               ", p = ", biostatUZH::formatPval(ma_meta_sere$pval.random),
               ", skewness of 0.00")
ESedgeadjust <- paste0(round(me95$estimate, 2), 
                       " (95\\% confidence interval from ", 
                       round(me95$CI[1], 2), " to ", 
                       round(me95$CI[2], 2),
                       ", p = ", biostatUZH::formatPval(me95$pval),
                       ", skewness of ", sprintf("%.2f",skew_pi(me95$CI, me95$estimate)))

## Overview COVID-19 Data (Table 1)
## -----------------------------------------------------------------------------
datHPH2025$steroids <- c("2 / 7", "69 / 128", "95 / 324", "11 / 75", "6 / 15", 
                         "26 / 105", "13 / 24")
datHPH2025$nosteroids <- c("2 / 12", "76 / 128", "283 / 683", "20 / 73",
                           "2 / 14", "29 / 92", "13 / 23")
datHPH2025$lower <- c(0.21, 0.49, 0.44, 0.20, 0.65, 0.38, 0.29)
datHPH2025$upper <- c(18.69, 1.31, 0.78, 1.04, 24.66, 1.33, 2.87)
datHPH2025$name <- c("DEXA-COVID 19", "CoDEX", "RECOVERY", "CAPE COVID", "COVID STEROID",
                     "REMAP-CAP", "Steroids-SARI")
sumcovid <- datHPH2025 |>
  mutate(
    `OR` = sprintf("%.2f", exp(logOR)),
    `95% CI` = sprintf("%.2f to %.2f", lower, upper)
  ) |>
  select(Study = name, Steroids = steroids, `No steroids` = nosteroids, OR, `95% CI`) |>
  xtable(
    caption = "Summary of seven randomized controlled trials investigating the 
    association between corticosteroids
    and mortality in hospitalized COVID-19
    patients \\citep{who2020corticosteroids}. Reported effect estimates are odds ratios.",
    label = "tab:coviddata",
    align = c("l", "l", "r", "r", "r", "r"),
    digits = c(0, 0, 0, 2, 2, 0)
  )

# add footnote
addtorow <- list()
addtorow$pos <- list(-1, nrow(datHPH2025))
addtorow$command <- c(
  "\\multicolumn{1}{l}{} & \\multicolumn{2}{c}{Deaths / Patients} & \\multicolumn{2}{c}{} \\\\ \n",
  "\\multicolumn{5}{l}{\\footnotesize CI = confidence interval, OR = odds ratio.} \\\\ \n"
)

# print table
print(sumcovid,
      include.rownames = FALSE,
      caption.placement = "top",
      hline.after = NULL, 
      booktabs = TRUE,
      add.to.row = addtorow)

## COVID-19 Example: Forest plot
## -----------------------------------------------------------------------------
# CD-Edgington confidence distribution
cdef <- remaeffect(escov, secov, "GAQ")$fcd
seqedge <- seq(me95$CI[1] - 0.4, me95$CI[2] + 0.4, length.out = 400)
edgecd <- data.frame(mu = seqedge, cd = cdef(seqedge))
edgecd$scaled <- edgecd$cd / max(edgecd$cd) * 0.75 # scale for plotting

# HKSJ confidence distribution
seqhksj <- seq(ma_meta_sere$lower.random - 0.4,
               ma_meta_sere$upper.random + 0.4, 
               length.out = 400)
hksj <- data.frame(
  mu = seqhksj,
  cd = dt((seqhksj - ma_meta_sere$TE.random) / ma_meta_sere$seTE.random,
          df = length(escov) - 1) / ma_meta_sere$seTE.random)
hksj$scaled <- hksj$cd / max(hksj$cd) * 0.75 # scale for plotting

# Dataset for plotting
fordat <- datHPH2025 |>
  mutate(
    l = logOR - 1.96 * logSE,
    u = logOR + 1.96 * logSE,
    t = "study"
  ) |>
  bind_rows(
    data.frame(
      name = c("CD-Edgington", "Hartung-Knapp-\nSidik-Jonkman"),
      logOR = c(me95$estimate, ma_meta_sere$TE.random),
      l = c(me95$CI[1], ma_meta_sere$lower.random),
      u = c(me95$CI[2], ma_meta_sere$upper.random),
      t = c("e1", "e2") 
    )
  ) |>
  mutate(
    name = factor(name, levels = c("Hartung-Knapp-\nSidik-Jonkman", "CD-Edgington", 
                                   rev(datHPH2025$name))),
    ynum = ifelse(as.numeric(name) > 2, as.numeric(name) / 2 + 1.5, as.numeric(name))
  )

# Plot
ggplot(fordat, aes(x = logOR, y = ynum, color = t)) +
  geom_errorbarh(aes(xmin = l, xmax = u), width = 0.2) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_line(
    data = edgecd, 
    aes(x = mu, y = 2 + scaled),
    color = "#F1C06E", 
    linewidth = 0.8
  ) +
  geom_line(
    data = hksj,
    aes(x = mu, y = 1 + scaled),
    color = "#90A4AE",
    linewidth = 0.8
  ) +
  scale_color_manual(values = c("#F1C06E", "#90A4AE", "black")) +
  scale_y_continuous(breaks = fordat$ynum, labels = fordat$name) +
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

## Covid-19 example: Computation of predictive distributions
##------------------------------------------------------------------------------
set.seed(021098)
# PCD-Fixed
pdsere1 <- PredDist(escov, secov, method = "FixedTau2")
pdsere1_pi <- pdsere1$PI

# PCD-simplified
pdsere2 <- PredDist(escov, secov, method = "SimplifiedCD")
pdsere2_pi <- pdsere2$PI

# PCD-full
pdsere3 <- PredDist(escov, secov, method = "FullCD")
pdsere3_pi <- pdsere3$PI

# PD: HTS
pdhts <- rt(1e7, df = length(escov) - 2) * # Student's t
  sqrt(ma_meta_sere$seTE.random^2 + ma_meta_sere$tau2) + # scale
  ma_meta_sere$TE.random # shift

## Confidence of theta_new >= 0 (for 4 displayed predictive distributions)
## -----------------------------------------------------------------------------
gr0 <- round(c(mean(pdhts >= 0), 
               mean(pdsere1$samples[, "theta_new"] >= 0),
               mean(pdsere2$samples[, "theta_new"] >= 0),
               mean(pdsere3$samples[, "theta_new"] >= 0)),
             3)

## Fisher skewness of predictive distributions distributions
## -----------------------------------------------------------------------------
skdis <- round(c(0, # HTS symmetric
                 EnvStats::skewness(pdsere1$samples[, "theta_new"]),
                 EnvStats::skewness(pdsere2$samples[, "theta_new"]),
                 EnvStats::skewness(pdsere3$samples[, "theta_new"])), 
               3)

## COVID-19 example: Display predictive distributions
## -----------------------------------------------------------------------------
# Combine samples from predictive distribution
sere_PDS <- data.frame( 
  hts = NA_real_, # display exact density rather than samples
  fixed = pdsere1$samples[, "theta_new"],
  simple = pdsere2$samples[, "theta_new"],
  full = pdsere3$samples[, "theta_new"]) |>
  pivot_longer(c(hts, fixed, simple, full), values_to = "s", names_to = "m") |>
  mutate(m = factor(m, levels = c("hts", "fixed", "simple", "full")))

# Quantile 95% and 99% prediction intervals
telePI <- data.frame(
  med = rep(c(median(pdhts), 
              median(pdsere1$samples[, "theta_new"]), 
              median(pdsere2$samples[, "theta_new"]), 
              median(pdsere3$samples[, "theta_new"])), 2),
  
  l = c(quantile(pdhts, 0.025), pdsere1$PI[1], pdsere2$PI[1], pdsere3$PI[1],
        quantile(pdhts, 0.005), quantile(pdsere1$samples[, "theta_new"], 0.005),
        quantile(pdsere2$samples[, "theta_new"], 0.005),
        quantile(pdsere3$samples[, "theta_new"], 0.005)),
  
  u = c(quantile(pdhts, 0.975), pdsere1$PI[2], pdsere2$PI[2], pdsere3$PI[2],
        quantile(pdhts, 0.995), quantile(pdsere1$samples[, "theta_new"], 0.995),
        quantile(pdsere2$samples[, "theta_new"], 0.995),
        quantile(pdsere3$samples[, "theta_new"], 0.995)),
  
  lev = rep(c("95%", "99%"), each = 4),
  
  m = factor(rep(c("hts", "fixed", "simple", "full"), 2),
             levels = c("hts", "fixed", "simple", "full"),
             labels = c("Higgins-Thompson-Spiegelhalter",
                        "PCD-fixed", "PCD-simplified", "PCD-full")),
  
  ypos = rep(c(2, 2.9, 1.9, 1.9), 2)
)

# 95% and 99% HPCD intervals for each distribution
hpd_95_pdhts   <- quantile(pdhts, c(0.025, 0.975)) # HTS distr. symmetric
hpd_95_pdsere1 <- HPDinterval(as.mcmc(pdsere1$samples[, "theta_new"]), prob = 0.95)
hpd_95_pdsere2 <- HPDinterval(as.mcmc(pdsere2$samples[, "theta_new"]), prob = 0.95)
hpd_95_pdsere3 <- HPDinterval(as.mcmc(pdsere3$samples[, "theta_new"]), prob = 0.95)

hpd_99_pdhts   <- quantile(pdhts, c(0.005, 0.995)) # HTS distr. symmetric
hpd_99_pdsere1 <- HPDinterval(as.mcmc(pdsere1$samples[, "theta_new"]), prob = 0.99)
hpd_99_pdsere2 <- HPDinterval(as.mcmc(pdsere2$samples[, "theta_new"]), prob = 0.99)
hpd_99_pdsere3 <- HPDinterval(as.mcmc(pdsere3$samples[, "theta_new"]), prob = 0.99)

telePI$l2 <- c(hpd_95_pdhts[1], hpd_95_pdsere1[1], hpd_95_pdsere2[1], hpd_95_pdsere3[1],
               hpd_99_pdhts[1], hpd_99_pdsere1[1], hpd_99_pdsere2[1], hpd_99_pdsere3[1])
telePI$u2 <- c(hpd_95_pdhts[2], hpd_95_pdsere1[2], hpd_95_pdsere2[2], hpd_95_pdsere3[2],
               hpd_99_pdhts[2], hpd_99_pdsere1[2], hpd_99_pdsere2[2], hpd_99_pdsere3[2])

# Determine break points histogram
brks <- seq(range(sere_PDS$s, na.rm = TRUE)[1],
            range(sere_PDS$s, na.rm = TRUE)[2],
            length.out = 500 + 1)
mids <- (brks[-length(brks)] + brks[-1]) / 2
bin_w <- brks[2] - brks[1]

# Density of HTS predictive distribution
df_hts <- length(escov) - 2
theta_hat <- ma_meta_sere$TE.random
se_theta <- sqrt(ma_meta_sere$seTE.random^2 + ma_meta_sere$tau2)
theta_grid <- seq(-2.5, 2.5, length.out = 1000)
htsdens <- data.frame(
  x = theta_grid,
  f = dt((theta_grid - theta_hat)/se_theta, df = df_hts) / se_theta,
  m = factor("hts", 
             levels = c("hts", "fixed", "simple", "full"),
             labels = c("Higgins-Thompson-Spiegelhalter",
                        "PCD-fixed", "PCD-simplified", "PCD-full")),
  region = ifelse(theta_grid < 0, "n", "p")
)

# Plot
sere_PDS |>
  mutate(
    m = factor(m, levels = c("hts", "fixed", "simple", "full"),
               labels = c("Higgins-Thompson-Spiegelhalter",
                          "PCD-fixed", "PCD-simplified", "PCD-full")),
    bin_idx = as.integer(cut(s, breaks = brks, include.lowest = TRUE, right = FALSE))) |>
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
  ungroup() |>
  ggplot(aes(x = bin_mid, y = density, fill = region, colour = region)) +
  facet_wrap(~ m, nrow = 2) +
  geom_col(width = bin_w, show.legend = FALSE) +
  geom_line(data = htsdens, inherit.aes = FALSE, aes(x = x, y = f, color = region)) +
  geom_area(data = htsdens, 
            aes(x = x, y = f, fill = region), 
            alpha = 1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("n" = "gray70", "p" = "#D98B8B"), guide = "none") +
  scale_color_manual(values = c("n" = "gray70", "p" = "#D98B8B"), guide = "none") +
  ggnewscale::new_scale_color() +
  geom_errorbarh(data = telePI |> mutate(type = "HCDP"),
                 aes(y = ypos, xmin = l, xmax = u,
                     alpha = lev, linewidth = lev, color = type),
                 inherit.aes = FALSE, height = 0) +
  geom_errorbarh(data = telePI |> mutate(type = "Equi-tailed"),
                 aes(y = ypos - 0.15, xmin = l2, xmax = u2,
                     alpha = lev, linewidth = lev, color = type),
                 inherit.aes = FALSE, height = 0) +
  geom_point(data = telePI,
             aes(x = med, y = ypos),
             inherit.aes = FALSE,
             size = 2, color ="#B88DCB") +
  geom_point(data = telePI,
             aes(x = med, y = ypos - 0.15),
             inherit.aes = FALSE,
             size = 2, color = "#6C91BF") +
  geom_segment(
    inherit.aes = FALSE,
    data = data.frame(
      m = factor(c("hts", "fixed", "simple", "full"),
                 levels = c("hts", "fixed", "simple", "full"),
                 labels = c("Higgins-Thompson-Spiegelhalter",
                            "PCD-fixed", "PCD-simplified", "PCD-full")),
      x = 1.5, y = 1, xend = c(0.15, 0.15, 0.2, 0.2), 
      yend = c(0.1, 0.1, 0.1, 0.1)
    ), aes(x = x, y = y, xend = xend, yend = yend)
  ) +
  geom_label(data = data.frame(
    m = factor(c("hts", "fixed", "simple", "full"),
                 levels = c("hts", "fixed", "simple", "full"),
                 labels = c("Higgins-Thompson-Spiegelhalter",
                            "PCD-fixed", "PCD-simplified", "PCD-full")),
    gr0 = gr0,
    x = 1.5,
    y = 1
  ), inherit.aes = FALSE, aes(x = x, y = y, label = gr0)) +
  scale_y_continuous(expand = c(0.0, 0), limits = c(0.0, 3.1)) +
  scale_alpha_manual(values = c(1, 0.75), guide = "none") +
  coord_cartesian(xlim = c(-2.5, 2.5)) +
  scale_linewidth_manual(values = c(1.5, 0.75)) +
  scale_color_manual(values = c("#6C91BF", "#B88DCB")) +
  labs(x = lab_tn, y = lab_ftn, linewidth = "Prediction interval", color = "") +
  theme_dk() +
  theme(strip.background = element_rect(fill = "gray99", color = "black"),
        legend.title = element_text(size = 8, face = "plain"))

## Covid Example: Prediction intervals (Table)
##------------------------------------------------------------------------------
# Function to find confidence of theta_new >= 0 for meta-methods
find_conf <- function(mm) {
  x <- optimize(function(lpi) {
    varest <- "PM"
    if (mm %in% c("KR", "KR-PR")) varest <- "REML"
    
    me <- metagen(escov, secov, random = T,
                  method.predict = mm,
                  method.tau = varest,
                  level.predict = lpi)
    abs(me$upper.predict)
  }, lower = 0, upper = 1)$minimum 
  
  (1 - x) / 2
}
# Find confidences for meta-methods
meta_pmethods <- c("V", "HTS", "HK", "HK-PR", "KR", "KR-PR", "NNF", "S")
gr0all <- unlist(lapply(meta_pmethods, find_conf))
gr0all[2] <- gr0[1] # exact for HTS by sampling (wo/ num. optim. error)

# 95% prediction intervals for all meta-methods
meta_combos <- expand.grid(method = meta_pmethods, level = 0.95)
sere_metapis <- as.data.frame(
  do.call(rbind, lapply(1:nrow(meta_combos), function(i) {
    method <- meta_combos$method[i]
    level  <- meta_combos$level[i]
    
    varest <- "PM"
    if (method %in% c("KR", "KR-PR")) varest <- "REML" # KR-PR requires REML estimation

    m <- metagen(
      TE = escov,
      seTE = secov,
      method.tau = varest,
      B = 100000,
      method.predict = method,
      level.predict = level
    )

    data.frame(
      method = method,
      level = paste0(level * 100, "%"),
      estimate = m$TE.random,
      pi.lower = m$lower.predict,
      pi.upper = m$upper.predict
    )
  }))
)

# Proposed methods
mp_pis <- data.frame(
  method = c("FixedTau2", "SimplifiedCD", "FullCD"),
  level = rep("95%", 3),
  estimate = telePI$med[2:4],
  pi.lower = c(pdsere1_pi[1], pdsere2_pi[1], pdsere3_pi[1]),
  pi.upper = c(pdsere1_pi[2], pdsere2_pi[2], pdsere3_pi[2])
)

# Bind prediction intervals together 
sere_pi <- sere_metapis |>
  bind_rows(mp_pis) |>
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
  sere_pi |>
    bind_cols(data.frame(gr = c(gr0all, (gr0)[2:4])))|>
    rowwise() |>
    mutate(
      sk = skew_pi(c(pi.lower, pi.upper), estimate),
      w = ifelse(pi.upper - pi.lower > 10^6, "$>10^6$", as.character(sprintf("%.2f", pi.upper - pi.lower))),
      method = factor(method,
      levels = c("PCD-full", "PCD-simplified", "PCD-fixed",
                 "Bootstrap", "HTS-KR-PR", "HTS-KR", "HTS-HK-PR",
                 "HTS-HK", "HTS-Veroniki", "HTS", "Skipka")),
      sk = case_when(
        round(sk, 4) == 0 ~ 0, 
        TRUE ~ sk
      ),
      estimate = ifelse(is.na(pi.lower), "\\textit{NC}", as.character(sprintf("%.2f", estimate))),
      gr = ifelse(is.na(pi.lower), " ", as.character(sprintf("%.3f", gr ))),
      pi = case_when(
        is.na(pi.lower) | is.na(pi.upper) ~ " ",
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
  align = c("l", "l", "r", "r", "r", "r", "r"),
  digits = c(0,0,2,0,2,3,3),
  caption = "95\\% prediction intervals from \\textsf{R} \\textsf{meta} package and Edgington’s PCD distributions, based on seven randomized controlled trials investigating the association between corticosteroids and mortality in hospitalized COVID-19 patients \\citep{who2020corticosteroids}. Medians of the predictive distributions and estimated confidence probabilities that a future effect is zero or larger are also shown. Note that the HTS-KR method requires the REML heterogeneity estimator, which is used for this approach.")

# add footnote
addtorow <- list()
addtorow$pos <- list(nrow(pixtabs))
addtorow$command <- c(" \n\\multicolumn{6}{l}{\\footnotesize HK = Hartung--Knapp, KR = Kenward--Roger, NC = non-convergence, PI = prediction interval, PR = Partlett--Riley.} \\\\ \n")

# print table
print(pixtabs,
      caption.placement = "top",
      include.rownames = FALSE,
      hline.after = NULL,
      booktabs = TRUE,
      add.to.row = addtorow,
      sanitize.text.function = identity)

## Simulation results: Load and preprocess results
##------------------------------------------------------------------------------
simresi <- list.files("../../results")
dt_v <- do.call(rbind, lapply(simresi, function(i) {
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

## Non-convergences of 95% prediction intervals
##------------------------------------------------------------------------------
dt <- dt_v[complete.cases(dt_v), ]
nc_nnf <- dt_v |>
  dplyr::select(starts_with(c("pi.cvr.", "ci.cvr.held.u"))) |>
  summarise(
    hts = mean(is.na(pi.cvr.hts)),
    nnf = mean(is.na(pi.cvr.nnf)),
    fix = mean(is.na(pi.cvr.fix)),
    simple = mean(is.na(pi.cvr.simple)),
    full = mean(is.na(pi.cvr.full)),
  ) |>
  pull(nnf)

## ReSim: Coverage of 95% Prediction Intervals
##------------------------------------------------------------------------------
# hline data
hlines <- expand.grid(
  k_large = paste0("Large~studies==", 0:2),
  I2 = "I^2==0*'%'"
)
hlines2 <- expand.grid(
  k_large = paste0("Large~studies==", 0:2),
  I2 = c("I^2==30*'%'", "I^2==60*'%'", "I^2==90*'%'")
)

simplot(dt, "N", "pi.cvr.", nam = 1) +
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

## ReSim: Skewness agreement of PIs with effect estimates
## -----------------------------------------------------------------------------
plotcor(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.hes, nam = 1) +
  labs(y = "Pearson correlation \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)

## ReSim: Width of 95% Prediction Intervals
##------------------------------------------------------------------------------
# Normal effects
simplot(dt = dt, distr = "N", stw_var = "pi.w.", nam = 1) +
  labs(y = "Width of 95% prediction intervals \u00B1 MCSE") +
  scale_y_log10()
# Warning message:
# Removed 60 rows containing missing values or values outside the scale range (`geom_point()`)
# This is due to the HTS interval, which is symmetric and has always zero skewness, i.e.
# NA, (60 conditions = 4I2 * 3largestudies * 5studies conditions)

################################################################################
################################################################################

## -----------------------------------------------------------------------------
## Code for Supplementary Material of Manuscript:
##  "Edgington's Method for Random-Effects Meta-Analysis Part I: Estimation"
## by
##   David Kronthaler and Leonhard Held
## -----------------------------------------------------------------------------

## Supplementary Table 1: Non-convergences parametric bootstrap PI
## -----------------------------------------------------------------------------
dtncnnf <- dt_v[!complete.cases(dt_v), ]
dtncnnf |>
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
    "Higgins' $I^2$ (\\%)" = "I2",
    "Distribution" = "dist"
  ) |> 
  xtable(label = "tab:ncnnf",
         align = c("r", "r", "r", "r", "l", "r"),
         caption = "Scenarios of the simulation study under which the 95\\% parametric bootstrap prediction interval did not converge at least once.") |>
  print(include.rownames = FALSE,
        sanitize.colnames.function = identity,
        caption.placement = "top",
        hline.after = NULL,
        booktabs = TRUE)

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

## ReSim: Coverage of 95% Prediction Intervals
##------------------------------------------------------------------------------
# hline data
hlines <- expand.grid(
  k_large = paste0("Large~studies==", 0:2),
  I2 = "I^2==0*'%'"
)
hlines2 <- expand.grid(
  k_large = paste0("Large~studies==", 0:2),
  I2 = c("I^2==30*'%'", "I^2==60*'%'", "I^2==90*'%'")
)
# Skew-normal effects
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

## ReSim: 95% PI coverage histograms for (Higgins I2 = 0%)
## -----------------------------------------------------------------------------
pcvr(0, NULL) +     
  labs(y = "Density", x = "Coverage probability") 

## ReSim: 95% PI coverage histograms for (Normal effects, Higgins I2 = 30%)
## -----------------------------------------------------------------------------
pcvr(30, "N")         

## ReSim: 95% PI coverage histograms for (Skew-normal effects, Higgins I2 = 30%)
## -----------------------------------------------------------------------------
pcvr(30, "LSN")         

## ReSim: 95% PI coverage histograms for (Normal effects, Higgins I2 = 60%)
## -----------------------------------------------------------------------------
pcvr(60, "N")         

## ReSim: 95% PI coverage histograms for (Skew-normal effects, Higgins I2 = 60%)
## -----------------------------------------------------------------------------
pcvr(60, "LSN")         

## ReSim: 95% PI coverage histograms for (Normal effects, Higgins I2 = 90%)
## -----------------------------------------------------------------------------
pcvr(90, "N")          

## ReSim: 95% PI coverage histograms for (Skew-normal effects, Higgins I2 = 90%)
## -----------------------------------------------------------------------------
pcvr(90, "LSN")          

## ReSim: Skewness agreement of PIs with effect estimates
## -----------------------------------------------------------------------------
# Skew-normal effects
plotcor(dt = dt, distr = "LSN", stw_var = "pi.sk.", var = sk.hes, nam = 1) +
  labs(y = "Pearson correlation \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning is due to HTS PI, see comment above

## ReSim: Skewness sign agreement of PIs with effect estimates (Cohens Kappa)
## -----------------------------------------------------------------------------
# Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.hes, nam = 1) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning is due to HTS PI, see comment above

## ReSim: Skewness sign agreement of PIs with effect estimates (Cohens Kappa)
## -----------------------------------------------------------------------------
# Skew-Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.hes, nam = 1) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning is due to HTS PI, see comment above

## ReSim: Skewness agreement of prediction intervals with true effects
## -----------------------------------------------------------------------------
# Normal effects
plotcor(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.es, nam = 1) +
  labs(y = "Pearson correlation \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) 
# Warning is due to HTS PI, see comment above

## ReSim: Skewness agreement of prediction intervals with true effects
## -----------------------------------------------------------------------------
# Skew-normal effects
plotcor(dt = dt, distr = "LSN", stw_var = "pi.sk.", var = sk.es, nam = 1) +
  labs(y = "Pearson correlation \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) 
# Warning is due to HTS PI, see comment above

## ReSim: Skewness sign agreement of PIs with true effects (Cohens Kappa)
## -----------------------------------------------------------------------------
# Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.es, nam = 1) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning is due to HTS PI, see comment above

## ReSim: Skewness sign agreement of PIs with true effects (Cohens Kappa)
## -----------------------------------------------------------------------------
# Skew-Normal effects
plotkappa(dt = dt, distr = "N", stw_var = "pi.sk.", var = sk.es, nam = 1) +
  labs(y = "Cohens kappa \u00B1 MCSE") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25)
# Warning is due to HTS PI, see comment above

## ReSim: Width of 95% Prediction Intervals
##------------------------------------------------------------------------------
# Skew-Normal effects
simplot(dt = dt, distr = "LSN", stw_var = "pi.w.", nam = 1) +
  labs(y = "Width of 95% prediction intervals \u00B1 MCSE") +
  scale_y_log10()

## ReSim: CRPS
##------------------------------------------------------------------------------
# Normal effects
simplot(dt = dt, distr = "N", stw_var = "crps.", nam = 2) +
  labs(y = "CRPS \u00B1 MCSE")

## ReSim: CRPS
##------------------------------------------------------------------------------
# Skew-Normal effects
simplot(dt = dt, distr = "LSN", stw_var = "crps.", nam = 2) +
  labs(y = "CRPS \u00B1 MCSE")

## ReSim: Computation time
##------------------------------------------------------------------------------
# Normal effects
simplot(dt, "N", "t.", 5) +
  labs(y = "Computation time (seconds) \u00B1 MCSE") +
  scale_y_continuous(trans = "log2") 

## ReSim: Computation time
##------------------------------------------------------------------------------
# Skew-normal effects
simplot(dt, "LSN", "t.", 5) +
  labs(y = "Computation time (seconds) \u00B1 MCSE") +
  scale_y_continuous(trans = "log2") 
