#Simulation code
#Zoe Rand and Eric Ward
#Last updated: 5 November 2026
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(MASS)
library(patchwork)

#simulation colors
pal <- c("#065BA1", "#677798")

#function to simulate data
simulate_data <- function(
  n_releases_before_treat = 300000,
  n_releases_after_treat = 1000000,
  n_releases_before_control = 300000,
  n_releases_after_control = 300000,
  baseline_recovery_rate = exp(-8.87), # survival x sampling rate for soos creek in Washington
  B_treatment_control = exp(0), # effect of treatment
  exp_B_after = exp(0), # how much does baseline rate change before:after for all
  phi = 4,
  n_years_before = 20,
  n_years_after = 10,
  seed = 1234,
  CV_control = 0.1, #interannual variability for control
  CV_treatment = 0.1 #interannual variabliity for treatment
) {
  set.seed(seed)

  dat <- expand_grid(
    treatment_control = c("Control", "Treatment"),
    brood_year = 1:(n_years_before + n_years_after)
  ) |>
    dplyr::mutate(
      period = ifelse(brood_year <= n_years_before, "Before", "After")
    )

  # dat_releases_before_treat<-rnorm(n_years_before, mean = n_releases_before_treat, sd = n_releases_before_treat * CV_treatment)
  # dat_releases_after_treat<-rnorm(n_years_after, mean = n_releases_after_treat, sd = n_releases_after_treat * CV_treatment)
  # dat_releases_before_control<-rnorm(n_years_before, mean = n_releases_before_control, sd = n_releases_before_control * CV_control)
  # dat_releases_after_control<-rnorm(n_years_after, mean = n_releases_after_control, sd = n_releases_after_control * CV_control)

  dat <- dat |>
    dplyr::mutate(
      period_f = relevel(factor(period), ref = "Before"),
      mean_rel = case_when(
        treatment_control == "Treatment" &
          period == "Before" ~ n_releases_before_treat,
        treatment_control == "Treatment" &
          period == "After" ~ n_releases_after_treat,
        treatment_control == "Control" &
          period == "Before" ~ n_releases_before_control,
        treatment_control == "Control" &
          period == "After" ~ n_releases_after_control
      ),
      cv = if_else(treatment_control == "Treatment", CV_treatment, CV_control),
      release_total = rnorm(dplyr::n(), mean = mean_rel, sd = mean_rel * cv),
      log_rel_total = log(release_total),
      # calculate the rate for each cell
      rate = baseline_recovery_rate *
        if_else(treatment_control == "Treatment", B_treatment_control, 1) *
        if_else(
          period == "After",
          exp_B_after,
          1
        ),
      mu = rate * release_total
    )

  # print(
  #   ggplot(dat) +
  #     geom_point(aes(x = brood_year, y = release_total)) +
  #     facet_grid(treatment_control ~ period)
  # )
  # simulate recovered values
  dat$n_recovered <- rnbinom(nrow(dat), mu = dat$mu, size = phi)
  offsets <- c(
    n_releases_before_control,
    n_releases_before_treat,
    n_releases_after_control,
    n_releases_after_treat
  )
  return(list(dat = dat, offsets = offsets))
}

#function to calculate BACI metrics
#Baci metrics based on equations from Chevalier et al. 2019

baci_calc <- function(m) {
  c(
    baci_contrast = (m[["IA"]] - m[["IB"]]) - (m[["CA"]] - m[["CB"]]),
    ci_contribution = abs(m[["IA"]] - m[["IB"]]) - abs(m[["CA"]] - m[["CB"]]),
    ci_divergence = abs(m[["IA"]] - m[["CA"]]) - abs(m[["IB"]] - m[["CB"]])
  )
}

#function fit model to simulated data
fit_sim_dat <- function(
  n_releases_before_treat = 300000,
  n_releases_after_treat = 1000000,
  n_releases_before_control = 300000,
  n_releases_after_control = 300000,
  baseline_recovery_rate = exp(-8.87),
  B_treatment_control = exp(0),
  exp_B_after = exp(0),
  phi = 4,
  n_years_before = 20,
  n_years_after = 10,
  seed = 1234,
  CV_control = 0.1, #interannual variability for control
  CV_treatment = 0.1 #interannual variabliity for treatment
) {
  out <- simulate_data(
    n_releases_before_treat,
    n_releases_after_treat,
    n_releases_before_control,
    n_releases_after_control,
    baseline_recovery_rate,
    B_treatment_control,
    exp_B_after,
    phi,
    n_years_before,
    n_years_after,
    seed,
    CV_control,
    CV_treatment
  )

  sim_dat <- out$dat

  fit <- glmmTMB(
    n_recovered ~ treatment_control *
      period_f +
      (1 | brood_year) +
      offset(log_rel_total),
    family = nbinom2,
    data = sim_dat
  )
  #print(summary(fit))

  #Calculate BACI metrics
  EMM <- emmeans(
    fit,
    ~ treatment_control * period_f,
    offset = 0,
    plotit = FALSE
  )
  EMM_df <- as_tibble(EMM)
  EMM_df$key <- paste0(
    #keys to match to the equations above
    ifelse(EMM_df$treatment_control == "Treatment", "I", "C"),
    ifelse(EMM_df$period_f == "After", "A", "B")
  )

  EMM_df <- EMM_df %>%
    mutate(
      release_mean = case_when(
        treatment_control == "Treatment" &
          period_f == "Before" ~ n_releases_before_treat,
        treatment_control == "Treatment" &
          period_f == "After" ~ n_releases_after_treat,
        treatment_control == "Control" &
          period_f == "Before" ~ n_releases_before_control,
        treatment_control == "Control" &
          period_f == "After" ~ n_releases_after_control
      ),
      log_releases = log(release_mean)
    )

  #brood year random effect
  sigma2 <- as.numeric(VarCorr(fit)$cond$brood_year)
  #point estimate
  est <- exp(EMM_df$emmean + EMM_df$log_releases + (sigma2 / 2))

  #simulate from multivariate normal to get uncertainty
  nsim <- 10000
  draws <- mvrnorm(nsim, mu = EMM_df$emmean, Sigma = vcov(EMM))
  #estimates for each draw
  draws2 <- apply(draws, 1, function(x) {
    exp(x + EMM_df$log_releases + (sigma2 / 2))
  })

  keys <- EMM_df$key
  #point estimates
  pt <- baci_calc(setNames(est, keys))
  #uncertainty
  sims <- apply(draws2, 2, function(x) baci_calc(setNames(x, keys)))

  res <- tibble(
    metric = names(pt),
    estimate = unname(pt),
    lcl = apply(sims, 1, quantile, 0.025),
    ucl = apply(sims, 1, quantile, 0.975)
  ) %>%
    mutate(
      pos = ifelse(estimate >= 0, 1, 0), # 1 if the estimate is positive
      signif = ifelse(lcl <= 0 & ucl >= 0, 0, 1)
    ) # 1 if the interval doesn't overlap 0

  #print(siginf)
  return(list(fit = fit, res = res))
}

#simulation settings
nsims <- 100 #number of simulations to run
release_options <- seq(1, 4, by = 0.1) #releases in the after-treatment are this * releases in the before
seeds <- sample(5000, nsims, replace = FALSE)

#plot result function
plot_sim_results <- function(sim_res, title) {
  plot1 <- ggplot(sim_res) +
    geom_line(aes(
      x = release_options,
      y = proportion,
      group = type,
      color = type
    )) +
    geom_point(
      aes(x = release_options, y = proportion, color = type, shape = type),
      alpha = 0.5
    ) +
    facet_wrap(~metric_group_f, nrow = 1) +
    scale_color_manual(
      values = pal,
      labels = c("Value is positive", "Confidence intervals \ndo not overlap 0")
    ) +
    scale_shape_manual(
      values = c(16, 4),
      labels = c("Value is positive", "Confidence intervals \ndo not overlap 0")
    ) +
    labs(
      x = "Treatment releases/control releases",
      y = "Proportion of simulations",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal() +
    ggtitle(title)
  return(print(plot1))
}


#Scenarios
#1) Recovery rate is the same before and after
list_of_fits <- list()
sim_fits <- list()
sim_signif_BACI <- matrix(
  NA,
  nrow = length(release_options),
  ncol = length(seeds)
)
sim_signif_contribution <- matrix(
  NA,
  nrow = length(release_options),
  ncol = length(seeds)
)
sim_signif_divergence <- matrix(
  NA,
  nrow = length(release_options),
  ncol = length(seeds)
)
pos_BACI <- matrix(NA, nrow = length(release_options), ncol = length(seeds))
pos_contribution <- matrix(
  NA,
  nrow = length(release_options),
  ncol = length(seeds)
)
pos_divergence <- matrix(
  NA,
  nrow = length(release_options),
  ncol = length(seeds)
)

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit <- fit_sim_dat(
      n_releases_after_treat = release_options[i] * 300000,
      seed = seeds[s],
      phi = 15
    )
    sim_fits[[i]] <- sim_fit$fit
    sim_signif_BACI[i, s] <- sim_fit$res$signif[
      sim_fit$res$metric == "baci_contrast"
    ]
    sim_signif_contribution[i, s] <- sim_fit$res$signif[
      sim_fit$res$metric == "ci_contribution"
    ]
    sim_signif_divergence[i, s] <- sim_fit$res$signif[
      sim_fit$res$metric == "ci_divergence"
    ]
    pos_BACI[i, s] <- sim_fit$res$pos[sim_fit$res$metric == "baci_contrast"]
    pos_contribution[i, s] <- sim_fit$res$pos[
      sim_fit$res$metric == "ci_contribution"
    ]
    pos_divergence[i, s] <- sim_fit$res$pos[
      sim_fit$res$metric == "ci_divergence"
    ]
  }
  list_of_fits[[s]] <- sim_fits
}

sim_signif_BACI_cumulative <- rowSums(sim_signif_BACI) / nsims
sim_signif_contribution_cumulative <- rowSums(sim_signif_contribution) / nsims
sim_signif_divergence_cumulative <- rowSums(sim_signif_divergence) / nsims
divergence_pos_cumulative <- rowSums(pos_divergence) / nsims
BACI_pos_cumulative <- rowSums(pos_BACI) / nsims
contribution_pos_cumulative <- rowSums(pos_contribution) / nsims
sim_res_table <- tibble(
  release_options = release_options,
  signif_BACI_contrast = sim_signif_BACI_cumulative,
  signif_CI_contribution = sim_signif_contribution_cumulative,
  signif_CI_divergence = sim_signif_divergence_cumulative,
  pos_BACI = BACI_pos_cumulative,
  pos_contribution = contribution_pos_cumulative,
  pos_divergence = divergence_pos_cumulative
) %>%
  pivot_longer(
    c(-release_options),
    names_to = "metric",
    values_to = "proportion"
  ) %>%
  mutate(
    metric_group = case_when(
      metric %in% c("signif_BACI_contrast", "pos_BACI") ~ "BACI contrast",
      metric %in%
        c("signif_CI_contribution", "pos_contribution") ~ "CI-contribution",
      metric %in% c("signif_CI_divergence", "pos_divergence") ~ "CI-divergence"
    )
  ) %>%
  mutate(
    type = ifelse(
      metric %in%
        c(
          "signif_BACI_contrast",
          "signif_CI_contribution",
          "signif_CI_divergence"
        ),
      "signif",
      "pos"
    )
  ) %>%
  mutate(
    metric_group_f = factor(
      metric_group,
      levels = c("BACI contrast", "CI-contribution", "CI-divergence")
    )
  )


p1 <- plot_sim_results(sim_res_table, "Baseline")


#2) there is a decline of exp(-1) in the recovery rate in the after period
list_of_fits2 <- list()
sim_fits2 <- list()
sim_signif2 <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit2 <- fit_sim_dat(
      n_releases_after_treat = release_options[i] * 300000,
      exp_B_after = exp(-1),
      seed = seeds[s],
      phi = 15
    )
    sim_fits2[[i]] <- sim_fit2$fit
    sim_signif2[i, s] <- sim_fit2$sig
  }
  list_of_fits2[[s]] <- sim_fits2
}

sim_signif_cumulative2 <- rowSums(sim_signif2) / nsims

p2 <- plot_sim_results(
  sim_signif_cumulative2,
  "Small decrease in recovery \n rate in the after period "
)


#3) decline in the recovery rate in the after is exp(-1.5)
list_of_fits3 <- list()
sim_fits3 <- list()
sim_signif3 <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit3 <- fit_sim_dat(
      n_releases_after_treat = release_options[i] * 300000,
      exp_B_after = exp(-1.5),
      seed = seeds[s],
      phi = 15
    )
    sim_fits3[[i]] <- sim_fit3$fit
    sim_signif3[i, s] <- sim_fit3$sig
  }
  list_of_fits3[[s]] <- sim_fits3
}

sim_signif_cumulative3 <- rowSums(sim_signif3) / nsims

p3 <- plot_sim_results(
  sim_signif_cumulative3,
  "Large decrease in recovery \n rate in the after period"
)


#4) increase in data overdispersion
#mid-level overdispersion (phi = 7)
list_of_fits4 <- list()
sim_fits4 <- list()
sim_signif4 <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit4 <- fit_sim_dat(
      n_releases_after_treat = release_options[i] * 300000,
      exp_B_after = exp(0),
      seed = seeds[s],
      phi = 7
    )
    sim_fits4[[i]] <- sim_fit4$fit
    sim_signif4[i, s] <- sim_fit4$sig
  }
  list_of_fits4[[s]] <- sim_fits4
}

sim_signif_cumulative4 <- rowSums(sim_signif4) / nsims

p4 <- plot_sim_results(
  sim_signif_cumulative4,
  "Middle overdispersion"
)

#higher overdispersion (phi = 3)
list_of_fits5 <- list()
sim_fits5 <- list()
sim_signif5 <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit5 <- fit_sim_dat(
      n_releases_after_treat = release_options[i] * 300000,
      exp_B_after = exp(0),
      seed = seeds[s],
      phi = 3
    )
    sim_fits5[[i]] <- sim_fit5$fit
    sim_signif5[i, s] <- sim_fit5$sig
  }
  list_of_fits5[[s]] <- sim_fits5
}

sim_signif_cumulative5 <- rowSums(sim_signif5) / nsims

p5 <- plot_sim_results(
  sim_signif_cumulative5,
  "High overdispersion"
)


#5) change in monitoring time in the after period
#5 years (the ones above are 10 years)
list_of_fits6 <- list()
sim_fits6 <- list()
sim_signif6 <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit6 <- fit_sim_dat(
      n_releases_after_treat = release_options[i] * 300000,
      exp_B_after = exp(0),
      seed = seeds[s],
      phi = 15,
      n_years_after = 5
    )
    sim_fits6[[i]] <- sim_fit6$fit
    sim_signif6[i, s] <- sim_fit6$sig
  }
  list_of_fits6[[s]] <- sim_fits6
}

sim_signif_cumulative6 <- rowSums(sim_signif6) / nsims

p6 <- plot_sim_results(
  sim_signif_cumulative6,
  "5 years in the \nafter period"
)


#2 years (the ones above are 10 years)
list_of_fits7 <- list()
sim_fits7 <- list()
sim_signif7 <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit7 <- fit_sim_dat(
      n_releases_after_treat = release_options[i] * 300000,
      exp_B_after = exp(0),
      seed = seeds[s],
      phi = 15,
      n_years_after = 2
    )
    sim_fits7[[i]] <- sim_fit7$fit
    sim_signif7[i, s] <- sim_fit7$sig
  }
  list_of_fits7[[s]] <- sim_fits7
}

sim_signif_cumulative7 <- rowSums(sim_signif7) / nsims

p7 <- plot_sim_results(
  sim_signif_cumulative7,
  "2 years in the \nafter period"
)


#6) Recovery rate is the same in before and after but lower than baseline
#small decrease
list_of_fits8 <- list()
sim_fits8 <- list()
sim_signif8 <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit8 <- fit_sim_dat(
      baseline_recovery_rate = exp(-9.3),
      n_releases_after_treat = release_options[i] * 300000,
      seed = seeds[s],
      phi = 15
    )
    sim_fits8[[i]] <- sim_fit8$fit
    sim_signif8[i, s] <- sim_fit8$sig
  }
  list_of_fits8[[s]] <- sim_fits8
}

sim_signif_cumulative8 <- rowSums(sim_signif8) / nsims

p8 <- plot_sim_results(sim_signif_cumulative, "Small decrease in recovery rate")

#large decrease
list_of_fits9 <- list()
sim_fits9 <- list()
sim_signif9 <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit9 <- fit_sim_dat(
      baseline_recovery_rate = exp(-9.6),
      n_releases_after_treat = release_options[i] * 300000,
      seed = seeds[s],
      phi = 15
    )
    sim_fits9[[i]] <- sim_fit9$fit
    sim_signif9[i, s] <- sim_fit9$sig
  }
  list_of_fits9[[s]] <- sim_fits9
}

sim_signif_cumulative9 <- rowSums(sim_signif9) / nsims

p9 <- plot_sim_results(sim_signif_cumulative, "Large decrease in recovery rate")


#plot simulations together
des <- "
AA##
BBCC
DDEE
FFGG
HHII"
pall <- p1 +
  p8 +
  p9 +
  p2 +
  p3 +
  p4 +
  p5 +
  p6 +
  p7 +
  plot_layout(design = des, axes = "collect")

ggsave(
  "figures/simualtion_results.png",
  pall,
  width = 6,
  height = 8,
  units = "in",
  dpi = 600
)
