#Simulation code
#Zoe Rand and Eric Ward
#Last updated: 2 September 2026
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(patchwork)

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
  seed = 1234
) {
  set.seed(seed)

  dat <- expand_grid(
    treatment_control = c("Control", "Treatment"),
    brood_year = 1:(n_years_before + n_years_after)
  ) |>
    dplyr::mutate(
      period = ifelse(brood_year <= n_years_before, "Before", "After")
    )

  dat <- dat |>
    dplyr::mutate(
      period_f = relevel(factor(period), ref = "Before"),
      release_total = case_when(
        treatment_control == "Treatment" &
          period == "Before" ~ n_releases_before_treat,
        treatment_control == "Treatment" &
          period == "After" ~ n_releases_after_treat,
        treatment_control == "Control" &
          period == "Before" ~ n_releases_before_control,
        treatment_control == "Control" &
          period == "After" ~ n_releases_after_control
      ),
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
  seed = 1234
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
    seed
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
  EMM <- emmip(
    fit,
    treatment_control ~ period_f,
    offset = log(out$offset),
    plotit = FALSE,
    CIs = TRUE
  )

  siginf <- ifelse(EMM$LCL[4] >= EMM$UCL[3], 1, 0)
  #print(siginf)
  return(list(fit = fit, sig = siginf))
}

#simulation settings
nsims <- 100 #number of simulations to run
release_options <- seq(1, 4, by = 0.1) #releases in the after-treatment are this * releases in the before
seeds <- sample(5000, nsims, replace = FALSE)

#plot result function
plot_sim_results <- function(sim_res, title) {
  plot1 <- ggplot() +
    geom_point(aes(x = release_options, y = sim_res)) +
    labs(
      x = "Treatment releases/control releases",
      y = "Proportion signficant"
    ) +
    theme_classic() +
    ggtitle(title)
  return(print(plot1))
}


#Scenarios
#1) Recovery rate is the same before and after
list_of_fits <- list()
sim_fits <- list()
sim_signif <- matrix(NA, nrow = length(release_options), ncol = length(seeds))

for (s in 1:length(seeds)) {
  for (i in 1:length(release_options)) {
    sim_fit <- fit_sim_dat(
      n_releases_after_treat = release_options[i] * 300000,
      seed = seeds[s],
      phi = 15
    )
    sim_fits[[i]] <- sim_fit$fit
    sim_signif[i, s] <- sim_fit$sig
  }
  list_of_fits[[s]] <- sim_fits
}

sim_signif_cumulative <- rowSums(sim_signif) / nsims

p1 <- plot_sim_results(sim_signif_cumulative, "Baseline")


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
