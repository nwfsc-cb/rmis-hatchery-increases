#RMIS GLMM model
#Zoe Rand

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(patchwork)


# Data wrangling ---------------------------------------------------------

#read in releases
rel <- read.csv("data/releases_2010_2026.csv") |>
  dplyr::select(-record_code) |>
  dplyr::rename(tag_code = tag_code_or_release_id)
rel$cwt_1st_mark_count[which(is.na(rel$cwt_1st_mark_count))] <- 0
rel$cwt_2nd_mark_count[which(is.na(rel$cwt_2nd_mark_count))] <- 0
rel$cwt_total <- rel$cwt_1st_mark_count + rel$cwt_2nd_mark_count

#read in total recoveries
rec <- read_csv("data/joined_data.csv")

#read in model data

dat <- readRDS("data/data_for_modeling.rds")

head(dat[[1]])

#naming by treatment hatchery so it's easier to work with
names(dat) <- c("Soos", "Quinault", "Samish", "Lewis", "Naselle", "Minter")

head(dat$Soos)

#note only focusing on Soos, Naselle, and Minter for now
#adding "after" designations to control hatcheries
dat$Soos$period[dat$Soos$release_year >= 2020] <- "After"
dat$Naselle$period[dat$Naselle$release_year >= 2020] <- "After"
dat$Minter$period[dat$Minter$release_year >= 2019] <- "After"


#the recoveries are by fish so need to summarise for the model
mod_dat_soos <- dat$Soos %>%
  group_by(
    brood_year,
    hatchery_location_name,
    fishery_region,
    treatment_control,
    period
  ) %>%
  summarise(n_recovered = floor(sum(number_cwt_estimated))) %>%
  ungroup()

mod_dat_naselle <- dat$Naselle %>%
  group_by(
    brood_year,
    hatchery_location_name,
    fishery_region,
    treatment_control,
    period
  ) %>%
  summarise(n_recovered = floor(sum(number_cwt_estimated))) %>%
  ungroup()

mod_dat_minter <- dat$Minter %>%
  group_by(
    brood_year,
    hatchery_location_name,
    fishery_region,
    treatment_control,
    period
  ) %>%
  summarise(n_recovered = floor(sum(number_cwt_estimated))) %>%
  ungroup()

#using release data to get total releases
rel_by_hatch <- rel %>%
  group_by(brood_year, hatchery_location_name) %>%
  summarise(release_total = sum(cwt_total)) %>%
  mutate(log_rel_total = log(release_total))

mod_dat_soos <- mod_dat_soos %>%
  left_join(rel_by_hatch, by = c("brood_year", "hatchery_location_name"))

mod_dat_naselle <- mod_dat_naselle %>%
  left_join(rel_by_hatch, by = c("brood_year", "hatchery_location_name"))

mod_dat_minter <- mod_dat_minter %>%
  left_join(rel_by_hatch, by = c("brood_year", "hatchery_location_name"))


#removing regions that don't have enough data
#Soos removing Alaska because only one control recovery in the after
mod_dat_soos <- mod_dat_soos %>% filter(fishery_region != "Alaska")

#Naselle removing Oregon because there's no control-after for that region
mod_dat_naselle <- mod_dat_naselle %>% filter(fishery_region != "Oregon")

#Minter
#removing California because there's only one data point
mod_dat_minter <- mod_dat_minter %>% filter(fishery_region != "California")
#removing oregon because there's no before data in the treatment hatchery
mod_dat_minter <- mod_dat_minter %>% filter(fishery_region != "Oregon")


# Plotting model data ----------------------------------------------------
#hatchery colors
hatch_col <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02')

#plotting function
plot_mod_dat <- function(dat, cutoff_yr) {
  plt2 <- ggplot(dat) +
    geom_line(aes(
      x = brood_year,
      y = n_recovered,
      group = hatchery_location_name_fac,
      color = hatchery_location_name_fac
    )) +
    geom_point(
      aes(x = brood_year, y = n_recovered, fill = log_rel_total),
      shape = 21,
      color = "transparent"
    ) +
    geom_vline(aes(xintercept = cutoff_yr), linetype = "dashed") +
    facet_grid(treatment_control ~ fishery_region) +
    labs(
      x = "Brood year",
      y = "Recoveries",
      fill = "Total releases (log)",
      color = "Hatchery"
    ) +
    scale_x_continuous(breaks = seq(2010, 2022, by = 4)) +
    scale_color_manual(values = hatch_col) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      strip.background.x = element_rect(fill = "white"),
      strip.background.y = element_blank()
    )
  print(plt2)
  return(plt2)
}

#plots
#making nice names for the hatcheries on the plot
mod_dat_soos$hatchery_location_name_fac <- factor(
  mod_dat_soos$hatchery_location_name,
  levels = c(
    "SOOS CREEK HATCHERY",
    "CLARKS CRK HATCHERY",
    "GORST CR REARING PND",
    "GROVERS CR HATCHERY",
    "ISSAQUAH HATCHERY",
    "VOIGHTS CR HATCHERY"
  ),
  labels = c(
    "Soos Creek",
    "Clarks Creek",
    "Gorst Creek",
    "Grovers Creek",
    "Issaquah",
    "Voights Creek"
  )
)

mod_dat_naselle$hatchery_location_name_fac <- factor(
  mod_dat_naselle$hatchery_location_name,
  levels = c("NASELLE HATCHERY", "FORKS CREEK HATCHERY", "NEMAH HATCHERY"),
  labels = c("Naselle", "Forks Creek", "Nemah")
)

mod_dat_minter$hatchery_location_name_fac <- factor(
  mod_dat_minter$hatchery_location_name,
  levels = c(
    "MINTER CR HATCHERY",
    "CLEAR CREEK HATCHERY",
    "KALAMA CR HATCHERY",
    "TUMWATER FALLS HATCHERY"
  ),
  labels = c("Minter Creek", "Clear Creek", "Kalama Creek", "Tumwater Falls")
)

soos_recs <- plot_mod_dat(mod_dat_soos, 2020) + ggtitle("Soos creek")
naselle_recs <- plot_mod_dat(mod_dat_naselle, 2020) + ggtitle("Naselle")
minter_recs <- plot_mod_dat(mod_dat_minter, 2019) + ggtitle("Minter creek")

#ggsave("figures/soos_recs.png", soos_recs, dpi = 600)
#ggsave("figures/naselle_recs.png", naselle_recs, dpi = 600)
#ggsave("figures/minter_recs.png", minter_recs, dpi = 600)

# Model ------------------------------------------------------------------

## Soos  ------------------------------------------------------------------

#making "before" the reference category for period
mod_dat_soos$period_f <- relevel(factor(mod_dat_soos$period), ref = "Before")

#making Washington the reference category for region (has the most recoveries)
mod_dat_soos$fishery_region_f <- relevel(
  factor(mod_dat_soos$fishery_region),
  ref = "Washington"
)

fit1 <- glmmTMB(
  n_recovered ~ treatment_control *
    period_f +
    fishery_region_f +
    (1 | brood_year) +
    offset(log_rel_total),
  family = nbinom2,
  data = mod_dat_soos
)
summary(fit1)

#trying with spatial effect
fit2 <- glmmTMB(
  n_recovered ~ treatment_control *
    period_f *
    fishery_region_f +
    (1 | brood_year) +
    offset(log_rel_total),
  family = nbinom2,
  data = mod_dat_soos
)
summary(fit2)
#fit 1 (without regional interaction) has a lower AIC--within 2 but without interaction is a simpler model

#diagnostics
fit1_simres <- simulateResiduals(fit1)
fit2_simres <- simulateResiduals(fit2)

plot(fit1_simres)
plot(fit2_simres)

plotResiduals(fit1_simres, form = model.frame(fit1)$treatment_control)
plotResiduals(fit1_simres, form = model.frame(fit1)$fishery_region_f)
plotResiduals(fit1_simres, form = model.frame(fit1)$period_f)
#plotResiduals(fit2_simres, form = model.frame(fit1)$brood_year)

#plotting results
plot_mod_results <- function(fit) {
  EMM_2 <- emmeans(
    fit,
    ~ treatment_control * period_f * fishery_region_f,
    offset = log(1000),
    type = "response"
  )
  print(plot(EMM_2))

  EMM_ip_2 <- emmip(
    fit,
    treatment_control ~ period_f | fishery_region_f,
    CIs = TRUE,
    plotit = FALSE
  )
  comparison_plot <- ggplot(EMM_ip_2) +
    geom_linerange(
      aes(x = xvar, ymin = LCL, ymax = UCL, color = tvar),
      linewidth = 2,
      alpha = 0.5
    ) +
    geom_line(aes(x = xvar, y = yvar, color = tvar, group = tvar)) +
    geom_point(aes(x = xvar, y = yvar, color = tvar)) +
    facet_wrap(~fishery_region_f) +
    labs(y = "Predicted recoveries (log scale)") +
    scale_color_manual(values = c("#1b9e77", "#7570b3")) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      strip.background = element_rect(fill = "white")
    )
  return(print(comparison_plot))
}

soos_plot <- plot_mod_results(fit1) + ggtitle("a) Soos creek")

## Naselle  ------------------------------------------------------------------

#making "before" the reference category for period
mod_dat_naselle$period_f <- relevel(
  factor(mod_dat_naselle$period),
  ref = "Before"
)

#making BC the reference category for region (has the most recoveries)
mod_dat_naselle$fishery_region_f <- relevel(
  factor(mod_dat_naselle$fishery_region),
  ref = "BC"
)


fit3 <- glmmTMB(
  n_recovered ~ treatment_control *
    period_f +
    fishery_region_f +
    (1 | brood_year) +
    offset(log_rel_total),
  family = nbinom2,
  data = mod_dat_naselle
)

summary(fit3)

fit4 <- glmmTMB(
  n_recovered ~ treatment_control *
    period_f *
    fishery_region_f +
    (1 | brood_year) +
    offset(log_rel_total),
  family = nbinom2,
  data = mod_dat_naselle
)
summary(fit4)
#model with interaction has lower AIC (fit 4)

#diagnostics
fit3_simres <- simulateResiduals(fit3)
fit4_simres <- simulateResiduals(fit4)


plot(fit3_simres)
plot(fit4_simres)

testQuantiles(fit3_simres)
testQuantiles(fit4_simres)

plotResiduals(fit4_simres, form = model.frame(fit4)$treatment_control)
plotResiduals(fit4_simres, form = model.frame(fit4)$fishery_region_f)
plotResiduals(fit4_simres, form = model.frame(fit4)$period_f)
#plotResiduals(fit4_simres, form = model.frame(fit4)$brood_year)

#some issues with the residuals currently

#plotting results
#prints two different ways of looking at the plots
naselle_plot <- plot_mod_results(fit4) + ggtitle("b) Naselle")

## Minter  ------------------------------------------------------------------
#making "before" the reference category for period
mod_dat_minter$period_f <- relevel(
  factor(mod_dat_minter$period),
  ref = "Before"
)

#making Washington the reference category for region (has the most recoveries)
mod_dat_minter$fishery_region_f <- relevel(
  factor(mod_dat_minter$fishery_region),
  ref = "Washington"
)


fit5 <- glmmTMB(
  n_recovered ~ treatment_control *
    period_f +
    fishery_region_f +
    (1 | brood_year) +
    offset(log_rel_total),
  family = nbinom2,
  data = mod_dat_minter
)

summary(fit5)


fit6 <- glmmTMB(
  n_recovered ~ treatment_control *
    period_f *
    fishery_region_f +
    (1 | brood_year) +
    offset(log_rel_total),
  family = nbinom2,
  data = mod_dat_minter
)

summary(fit6)
#model without interaction has lower AIC

#diagnostics
fit5_simres <- simulateResiduals(fit5)
fit6_simres <- simulateResiduals(fit6)


plot(fit5_simres)
plot(fit6_simres)

testQuantiles(fit5_simres)

plotResiduals(fit5_simres, form = model.frame(fit5)$treatment_control)
plotResiduals(fit5_simres, form = model.frame(fit5)$fishery_region_f)
plotResiduals(fit5_simres, form = model.frame(fit5)$period_f)


#some issues with the residuals currently

#plotting results
minter_plot <- plot_mod_results(fit5) + ggtitle("c) Minter creek")


#plotting all three together
plot_tog <- soos_plot /
  naselle_plot /
  minter_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

plot_tog

ggsave(
  "figures/initial_mod_results.png",
  plot_tog,
  dpi = 300,
  width = 6,
  height = 10,
  units = "in"
)
