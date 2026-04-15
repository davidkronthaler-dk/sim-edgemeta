## P-value function Wald test
##------------------------------------------------------------------------------
p_fun_wald <- function(h0_point, es, se, onesided = T){
  if (onesided) p <- 1-pnorm(es, h0_point, se)
  else if (!onesided) p <- 2 * apply(cbind(pnorm(es, h0_point, se), 
                                           1 - pnorm(es, h0_point, se)), 1, min)
  return(p)
}

## Function to calculate skewness of an interval
##------------------------------------------------------------------------------
skew_pi <- function(pi, point) (pi[2] + pi[1] - 2 * point)/(pi[2] - pi[1])

## Fishers weighted skewness coefficient
##------------------------------------------------------------------------------
fwskew <- function(es, se){
  wi <- 1/(se^2)
  mu_bar <- sum(es * wi)/sum(wi)
  num <- (sum(wi * (es -mu_bar)^3)) * sqrt(sum(wi))
  denum <- (sum(wi * (es - mu_bar)^2))^(3/2)
  return(num/denum)
}

## Mean with MCSE errorbars plot function for simulation results
## -----------------------------------------------------------------------------
simplot <- function(dt, distr, stw_var, nam, hts = TRUE) {
  
  # filter not displayed methods (for fixed-effect meta-analysis)
  dt <- dt |> dplyr::select(-bias.ivw.fix, - bias.held.u,
                            -sqe.ivw.fix, -sqe.held.u,
                            -ci.cvr.held.u, -ci.w.held.u)
  
  lev <- switch (nam,
    c("hts", "nnf", "fix", "simple", "full"),
    c("hts", "fix", "simple", "full"),
    c("ivw.random", "held.a", "samp"),
    c("r", "hk", "held.a", "samp"),
    c("h", "nnf", "fix", "simple", "full")
  )
  
  lab <- switch(nam,
                as.expression(c(
                  "Higgins-Thompson-Spiegelhalter",
                  "Bootstrap",
                  expression("PCD-fixed"),
                  "PCD-simplified",
                  "PCD-full"
                )),
                as.expression(c(
                  "Higgins-Thompson-Spiegelhalter",
                  expression("PCD-fixed"),
                  "PCD-simplified",
                  "PCD-full"
                )),
                as.expression(c("Random-effects", "Edgington (add. heterogeneity)",
                                "CD-Edgington")),
                as.expression(c("Random-effects", "Hartung-Knapp-Sidik-Jonkman",  
                                "Edgington (add. heterogeneity)", "CD-Edgington")),
                as.expression(c(
                  "Higgins-Thompson-Spiegelhalter",
                  "Bootstrap",
                  expression("PCD-fixed"),
                  "PCD-simplified",
                  "PCD-full"
                ))
  )
  
  if (!hts) dt <- dt |> dplyr::select(-pi.w.hts)
  
  cst <- c(
    hts = "#7C878EFF",
    nnf = "#84BD00FF",
    fix = "#FFCD00FF",
    simple = "#5C88DAFF",
    full = "#CC0C00FF",
    ivw.random = "#FFCD00FF",
    held.a = "#5C88DAFF",
    samp = "#CC0C00FF",
    r = "#7C878EFF",
    hk = "#FFCD00FF",
    h = "#7C878EFF"
  )
  cst <- cst[lev]
  
  dt |>
    filter(dist == distr) |>
    pivot_longer(starts_with(!!stw_var), values_to = "v", names_to = "m") |>
    mutate(
      m = factor(gsub(!!stw_var, "", m), levels = lev),
      I2 = factor(paste0("I^2==", I2, "*'%'"),
                  levels = c("I^2==0*'%'", "I^2==30*'%'",
                             "I^2==60*'%'", "I^2==90*'%'")),
      k = factor(k, levels = c(3,5,10,20,50)),
      k_large = factor(paste0("Large~studies==", k_large),
                       levels = paste0("Large~studies==", 0:2)),
    ) |>
    mutate(
      I2 = droplevels(I2)
    ) |>
    group_by(k, I2, k_large, m) |>
    summarise(
      me = mean(v, na.rm = T),
      se = sd(v, na.rm = T) / sqrt(sum(!is.na(v)) - 1),
      .groups = "drop"
    ) |>
    complete(k, I2, k_large, m, fill = list(me = NA, se = NA)) |> 
    ggplot(aes(x = k, y = me, group = m, color = m, fill = m)) +
    facet_grid(k_large ~ I2, drop = FALSE, labeller = label_parsed) +
    geom_point(position = position_dodge(0.4), size = 1.2, shape = 18) +
    geom_errorbar(aes(ymin = me - se, ymax = me + se), 
                  position = position_dodge(0.4), size = 0.5,
                  width = 0.75) +
    labs(x = "Number of studies") +
    scale_color_manual(
      values =  cst,
      labels = lab
    ) +
    scale_fill_manual(
      values = cst,
      labels = lab
    ) +
    theme_dk() +
    theme(
      axis.text = element_text(colour = "black", face = "italic",
                               family = "Times",
                               size = 8),
      axis.title = element_text(colour = "black",
                                size = 8),
      axis.ticks = element_line(colour = "black"),
      title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      strip.text = element_text(size = 8)
    )
}

## Errobarplot for correlation results from simulation study
##------------------------------------------------------------------------------
plotcor <- function(dt, distr, stw_var, var, nam) {
  
  # filter not displayed methods (fixed-effect meta-analysis)
  dt <- dt |> dplyr::select(-ci.sk.held.u)
  
  lev <- switch (nam ,
                 c("nnf", "fix", "simple", "full"),
                 c("held.a", "samp")
  )
  
  lab <- switch(nam,
                as.expression(c(
                  "Bootstrap",
                  expression("PCD-fixed"),
                  "PCD-simplified",
                  "PCD-full"
                )),
                as.expression(c("Edgington (add. heterogeneity)", 
                                "CD-Edgington"))
  )
  
  del <- if (deparse(substitute(var)) == "sk.es") TRUE else FALSE
  
  if (del) {
    dt <- dt |> filter(I2 != 0)
  }
  
  cst <- c(
    hts = "#7C878EFF",
    nnf = "#84BD00FF",
    fix = "#FFCD00FF",
    simple = "#5C88DAFF",
    full = "#CC0C00FF",
    ivw.random = "#FFCD00FF",
    held.a = "#5C88DAFF",
    samp = "#CC0C00FF",
    r = "#7C878EFF",
    hk = "#FFCD00FF",
    h = "#7C878EFF"
  )
  cst <- cst[lev]

  dt |>
    filter(dist == distr) |>
    pivot_longer(starts_with(!!stw_var), values_to = "v", names_to = "m") |>
    mutate(
      m = factor(gsub(!!stw_var, "", m), levels = lev), # labels = lab),
      I2 = factor(paste0("I^2==", I2, "*'%'"),
                  levels = c("I^2==0*'%'", "I^2==30*'%'",
                             "I^2==60*'%'", "I^2==90*'%'")),
      k = factor(k, levels = c(3,5,10,20,50)),
      k_large = factor(paste0("Large~studies==", k_large),
                       levels = paste0("Large~studies==", 0:2)),
    ) |>
    mutate(
      I2 = droplevels(I2)
    ) |>
    group_by(k, I2, k_large, m) |>
    summarise(
      cr = cor(v, {{var}}, use = "pairwise.complete.obs"),
      n_pairs = sum(!is.na(v) & !is.na({{var}})),
      se = sqrt((1 - cr^2) / (n_pairs - 2)),
      .groups = "drop"
    ) |>
    complete(k, I2, k_large, m, fill = list(me = NA, se = NA)) |> 
    ggplot(aes(x = k, y = cr, group = m, color = m, fill = m)) +
    facet_grid(k_large ~ I2, drop = FALSE, labeller = label_parsed) +
    geom_point(position = position_dodge(0.3), size = 1.2, shape = 18) +
    geom_errorbar(aes(ymin = cr - se, ymax = cr + se),
                  position = position_dodge(0.3), size = 0.5,
                  width = 0.75) +
    labs(x = "Number of studies") +
    scale_color_manual(
      values = cst,
      labels = lab,
      na.translate = FALSE
    ) +
    scale_fill_manual(
      values = cst,
      labels = lab,
      na.translate = FALSE
    ) +
    theme_dk() +
    theme(
      axis.text = element_text(colour = "black", face = "italic",
                               family = "Times",
                               size = 8),
      axis.title = element_text(colour = "black",
                                size = 8),
      axis.ticks = element_line(colour = "black"),
      title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      strip.text = element_text(size = 8)
    )
}

## Plot coverage probability distributions
## -----------------------------------------------------------------------------
pcvr <- function(i2f, distf) {
  
  lev <- c("hts", 
           "nnf", 
           "fix", 
           "simple", 
           "full")
  
  lab <- as.expression(c("HTS", 
                         "Bootstrap", 
                         expression("PCD-fixed"),
                         "PCD-simplified", 
                         "PCD-full"))
  
  cst <- c(hts = "#7C878EFF",
           nnf = "#84BD00FF", 
           fix = "#FFCD00FF",
           simple = "#5C88DAFF",
           full = "#CC0C00FF")
  
  I2_levels <- c("I^2==0*'%'", "I^2==30*'%'", "I^2==60*'%'", "I^2==90*'%'")
  k_levels <- paste0("k==", c(3, 5, 10, 20, 50))
  k_large_levels <- paste0("k[large]==", 0:2)
  
  vl <- ifelse(i2f == 0, 1, 0.95)
  
  if (!is.null(distf)) {
    dtx <- dt |>
      filter(I2 == i2f & dist == distf)
  } else {
    dtx <- dt |>
      filter(I2 == i2f)
  }
  
  dtx |>
    pivot_longer(cols = starts_with("pi.cvr."), values_to = "c", names_to = "m") |>
    mutate(
      m = factor(gsub("pi.cvr.", "", m), levels = lev),
      m_lab = factor(m, levels = lev, labels = lab), 
      I2 = factor("I^2==0*'%'", levels = I2_levels),
      k = factor(paste0("k==", k), levels = k_levels),
      k_large = factor(paste0("k[large]==", k_large), levels = k_large_levels)
    ) |>
    ggplot(aes(x = c, color = m)) + 
    facet_grid(m_lab ~ k + k_large, labeller = label_parsed) +
    geom_histogram(aes(y = after_stat(density)))  +
    scale_x_continuous(breaks = c(0, 0.5, vl), 
                       labels = c("0", "0.5", ifelse(vl == 1, "1", "0.95"))) +
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = cst, labels = lab) + 
    theme_dk() +
    theme(axis.text = element_text(colour = "black", face = "italic",
                                   family = "Times", size = 6),
      legend.position = "none",
      axis.title = element_text(colour = "black",
                                size = 8),
      axis.ticks = element_line(colour = "black"),
      title = element_text(size = 8, face = "bold"),
      strip.text = element_text(size = 6)
    ) +
    labs(y = "Density", x = "Coverage probability") +
    geom_vline(xintercept = vl, linetype = "dashed", size = 0.25)
}

## Cohens Kappa for sign agreement
## -----------------------------------------------------------------------------

plotkappa <- function(dt, distr, stw_var, var, nam) {
  
  lev <- switch (nam ,
                 c("nnf", "fix", "simple", "full"),
                 c("held.a", "samp")
  )
  
  lab <- switch(nam,
                as.expression(c(
                  "Bootstrap",
                  expression("PCD-fixed"),
                  "PCD-simplified",
                  "PCD-full"
                )),
                as.expression(c("Edgington (add. heterogeneity)", 
                                "CD-Edgington"))
  )
  
  del <- if (deparse(substitute(var)) == "sk.es") TRUE else FALSE
  
  if (del) {
    dt <- dt |> filter(I2 != 0)
  }
  
  cst <- c(
    hts = "#7C878EFF",
    nnf = "#84BD00FF",
    fix = "#FFCD00FF",
    simple = "#5C88DAFF",
    full = "#CC0C00FF",
    ivw.fix = "#7C878EFF",
    ivw.random = "#FFCD00FF",
    held.a = "#5C88DAFF",
    samp = "#CC0C00FF",
    r = "#7C878EFF",
    hk = "#FFCD00FF",
    h = "#7C878EFF"
  )
  cst <- cst[lev]
  
  dt |>
    filter(dist == distr) |>
    pivot_longer(starts_with(!!stw_var), values_to = "v", names_to = "m") |>
    mutate(
      m = factor(gsub(!!stw_var, "", m), levels = lev), # labels = lab),
      I2 = factor(paste0("I^2==", I2, "*'%'"),
                  levels = c("I^2==0*'%'", "I^2==30*'%'",
                             "I^2==60*'%'", "I^2==90*'%'")),
      k = factor(k, levels = c(3,5,10,20,50)),
      k_large = factor(paste0("Large~studies==", k_large),
                       levels = paste0("Large~studies==", 0:2)),
    ) |>
    mutate(
      I2 = droplevels(I2),
      sign_v = sign(v),
      sign_var = sign({{var}})
    ) |>
    group_by(k, I2, k_large, m) |>
    summarise(
      kappa = psych::cohen.kappa(cbind(sign_v, sign_var))$kappa,
      se = sqrt(psych::cohen.kappa(cbind(sign_v, sign_var))$var.kappa),
      .groups = "drop"
    ) |>
    complete(k, I2, k_large, m, fill = list(me = NA, se = NA)) |> 
    ggplot(aes(x = k, y = kappa, group = m, color = m, fill = m)) +
    facet_grid(k_large ~ I2, drop = FALSE, labeller = label_parsed) +
    geom_point(position = position_dodge(0.3), size = 1.2, shape = 18) +
    geom_errorbar(aes(ymin = kappa - se, ymax = kappa + se),
                  position = position_dodge(0.3), size = 0.5,
                  width = 0.75) +
    labs(x = "Number of studies") +
    scale_color_manual(
      values = cst,
      labels = lab,
      na.translate = FALSE
    ) +
    scale_fill_manual(
      values = cst,
      labels = lab,
      na.translate = FALSE
    ) +
    theme_dk() +
    theme(
      axis.text = element_text(colour = "black", face = "italic",
                               family = "Times",
                               size = 8),
      axis.title = element_text(colour = "black",
                                size = 8),
      axis.ticks = element_line(colour = "black"),
      title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      strip.text = element_text(size = 8)
    )
}
