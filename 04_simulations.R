#Simulation code
#Zoe Rand and Eric Ward
#Last updated: 5 November 2026
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(future.apply) # for parallelization
library(patchwork)


# Settings ---------------------------------------------------------------
nsims <- 1000 #number of simulated datasets per scenario and release level
release_options <- seq(1, 4, by = 0.1) #amount that releases increase in the after period from the before period

base_releases <- list(
  n_releases_before_treat = 300000,
  n_releases_before_control = 300000,
  n_releases_after_control = 300000
)

n_draws <- 10000 #draws from MVN for BACI metric uncertainty

set.seed(123)
seeds <- sample(5000, nsims, replace = FALSE)


#simulation colors
pal <- c("#065BA1", "#677798")

# Scenarios --------------------------------------------------------------
scenarios <- tribble(
  ~scenario             , ~title                                                      , ~args                                                                                 ,
  # Recovery rate same before/after
  "baseline"            , "a) Baseline"                                               , list(phi = 15)                                                                        ,
  # Higher releases in the control:
  "control_600k"        , "b) Higher releases in control"                             , list(phi = 15, n_releases_before_control = 600000, n_releases_after_control = 600000) ,
  # Recovery rate same before/after, but lower than baseline
  "small_rate_decrease" , "c) Small decrease in recovery rate"                        , list(baseline_recovery_rate = exp(-9.3), phi = 15)                                    ,
  "large_rate_decrease" , "d) Large decrease in recovery rate"                        , list(baseline_recovery_rate = exp(-9.6), phi = 15)                                    ,
  # Recovery rate declines in the after period
  "small_after_decline" , "e) Small decline in recovery \n rate in the after period " , list(exp_B_after = exp(-1), phi = 15)                                                 ,
  "large_after_decline" , "f) Large decline in recovery \n rate in the after period"  , list(exp_B_after = exp(-1.5), phi = 15)                                               ,
  # Overdispersion
  "mid_overdispersion"  , "g) Middle overdispersion"                                  , list(phi = 7)                                                                         ,
  "high_overdispersion" , "h) High overdispersion"                                    , list(phi = 3)                                                                         ,
  #More variability in control only
  "control_highCV"      , "i) High CV in control releases"                            , list(phi = 15, CV_control = 0.3)                                                      ,
  #more variability in both
  "highCV_all"          , "j) High CV in both control \nand treatment releases"       , list(phi = 15, CV_control = 0.3, CV_treatment = 0.3)                                  ,
  # Monitoring time in the after period (other scenarios use 10 years)
  "after_5yr"           , "k) 5 years in the after period"                            , list(phi = 15, n_years_after = 5)                                                     ,
  "after_2yr"           , "l) 2 years in the after period"                            , list(phi = 15, n_years_after = 2)
)


# Functions --------------------------------------------------------------
#function to simulate dataset
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

  return(list(dat = dat))
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
fit_sim_dat <- function(..., n_draws = n_draws) {
  out <- simulate_data(...)

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
  EMM_df <- as_tibble(EMM) %>%
    left_join(
      distinct(sim_dat, treatment_control, period_f, mean_rel),
      by = c("treatment_control", "period_f")
    )
  EMM_df$key <- paste0(
    #keys to match to the equations above
    ifelse(EMM_df$treatment_control == "Treatment", "I", "C"),
    ifelse(EMM_df$period_f == "After", "A", "B")
  )

  #brood year random effect
  sigma2 <- as.numeric(VarCorr(fit)$cond$brood_year)
  #point estimate
  est <- exp(EMM_df$emmean + log(EMM_df$mean_rel) + (sigma2 / 2))

  #simulate from multivariate normal to get uncertainty
  nsim <- 10000
  draws <- MASS::mvrnorm(nsim, mu = EMM_df$emmean, Sigma = vcov(EMM))
  #estimates for each draw
  draws2 <- apply(draws, 1, function(x) {
    exp(x + log(EMM_df$mean_rel) + (sigma2 / 2))
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
  return(res)
}

#function to run a simulation
#returns NULL if the fit fails so that a single bad fit doesn't stop the whole run
run_one <- function(args, release_mult, seed) {
  # scenario args override the defaults for each release number
  p <- modifyList(base_releases, args)
  p$n_releases_after_treat <- release_mult * p$n_releases_after_control
  p$seed <- seed
  p$n_draws <- n_draws

  #print(p)
  tryCatch(
    do.call(fit_sim_dat, p) |> dplyr::select(metric, pos, signif),
    error = function(e) NULL
  )
}

#plot result function
plot_sim_results <- function(sim_res, title) {
  plot1 <- ggplot(sim_res) +
    geom_line(aes(
      x = release_mult,
      y = proportion,
      group = type,
      color = type
    )) +
    geom_point(
      aes(x = release_mult, y = proportion, color = type, shape = type),
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
      x = "Treatment: After releases/Before releases",
      y = "Proportion of simulations",
      color = NULL,
      shape = NULL
    ) +
    theme_minimal() +
    ggtitle(title)
  return(print(plot1))
}


# Run scenarios ----------------------------------------------------------
scenario_grid <- expand_grid(
  scenarios,
  release_mult = release_options,
  seed = seeds
)
plan(multisession)
sim_results <- future_mapply(
  run_one,
  scenario_grid$args,
  scenario_grid$release_mult,
  scenario_grid$seed,
  future.seed = NULL,
  SIMPLIFY = TRUE
)

#combine results
sim_results_df <- t(sim_results) %>%
  as_tibble() %>%
  rowid_to_column(var = "Sim_ID") %>%
  unnest_longer(2:4) %>%
  dplyr::select(-"signif_id")
res <- scenario_grid %>%
  as_tibble() %>%
  dplyr::select(-args) %>%
  rowid_to_column(var = "Sim_ID") %>%
  left_join(sim_results_df, by = "Sim_ID")

#summarise results
res_df <- res %>%
  pivot_longer(c(signif, pos), names_to = "type", values_to = "value") %>%
  group_by(scenario, release_mult, metric, type) %>%
  summarise(proportion = sum(value) / nsims) %>%
  mutate(
    metric_group_f = factor(
      metric,
      levels = c("baci_contrast", "ci_contribution", "ci_divergence"),
      labels = c("BACI contrast", "CI-contribution", "CI-divergence")
    )
  )


# Plot results -----------------------------------------------------------
plts <- list()
for (s in 1:nrow(scenarios)) {
  plts[[s]] <- plot_sim_results(
    res_df[res_df$scenario == scenarios$scenario[s], ],
    scenarios$title[s]
  )
}

p_all <- wrap_plots(plts) +
  plot_layout(
    nrow = nrow(scenarios) / 2,
    axes = "collect",
    guides = "collect"
  ) &
  theme(legend.position = "bottom")
p_all

ggsave(
  "figures/simulation_results.pdf",
  p_all,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600
)
