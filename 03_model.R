#RMIS GLMM
#Zoe Rand
#last updated 2 September 2026

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(patchwork)
library(broom.mixed)
library(modelsummary) #for supplement tables


# Data wrangling ---------------------------------------------------------

#read in releases
rel <- read.csv("data/releases_2010_2026.csv") |>
  dplyr::select(-record_code) |>
  dplyr::rename(tag_code = tag_code_or_release_id)
rel$cwt_1st_mark_count[which(is.na(rel$cwt_1st_mark_count))] <- 0
rel$cwt_2nd_mark_count[which(is.na(rel$cwt_2nd_mark_count))] <- 0
rel$cwt_total <- rel$cwt_1st_mark_count + rel$cwt_2nd_mark_count

rel$non_cwt_1st_mark_count[which(is.na(rel$non_cwt_1st_mark_count))] <- 0
rel$non_cwt_2nd_mark_count[which(is.na(rel$non_cwt_2nd_mark_count))] <- 0
rel$total_releases <- rel$cwt_total +
  rel$non_cwt_1st_mark_count +
  rel$non_cwt_2nd_mark_count

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
dat$Soos$release_year <- as.numeric(substr(dat$Soos$first_release_date, 1, 4))
dat$Naselle$release_year <- as.numeric(substr(
  dat$Naselle$first_release_date,
  1,
  4
))
dat$Minter$release_year <- as.numeric(substr(
  dat$Minter$first_release_date,
  1,
  4
))
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
  summarise(release_total = sum(total_releases)) %>%
  mutate(log_rel_total = log(release_total))

mod_dat_soos <- mod_dat_soos %>%
  left_join(rel_by_hatch, by = c("brood_year", "hatchery_location_name"))

mod_dat_naselle <- mod_dat_naselle %>%
  left_join(rel_by_hatch, by = c("brood_year", "hatchery_location_name"))

mod_dat_minter <- mod_dat_minter %>%
  left_join(rel_by_hatch, by = c("brood_year", "hatchery_location_name"))

#removing hatcheries that aren't treatments or controls
mod_dat_soos <- mod_dat_soos %>%
  filter(hatchery_location_name != "CLARKS CRK HATCHERY")
mod_dat_naselle <- mod_dat_naselle %>%
  filter(hatchery_location_name != "FORKS CREEK HATCHERY")

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
plot_mod_dat <- function(dat, cutoff_yr, plot_title) {
  max_year <- max(dat$brood_year)
  print(max_year)
  plt2 <- ggplot(dat) +
    geom_line(aes(
      x = brood_year,
      y = n_recovered,
      group = hatchery_location_name_fac,
      color = hatchery_location_name_fac
    )) +
    geom_point(
      aes(x = brood_year, y = n_recovered, color = hatchery_location_name_fac)
      #shape = 21,
      #color = "transparent"
    ) +
    geom_vline(aes(xintercept = cutoff_yr), linetype = "dashed") +
    facet_grid(treatment_control ~ fishery_region) +
    labs(
      x = "Brood year",
      y = "CWT Recoveries",
      color = "Hatchery"
    ) +
    scale_x_continuous(breaks = seq(2010, max_year, by = 4)) +
    scale_color_manual(values = hatch_col) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(plot_title)

  release_dat <- dat %>%
    select(
      brood_year,
      release_total,
      hatchery_location_name_fac,
      treatment_control
    ) %>%
    distinct()
  plt3 <- ggplot(release_dat) +
    geom_line(aes(
      x = brood_year,
      y = release_total,
      group = hatchery_location_name_fac,
      color = hatchery_location_name_fac
    )) +
    geom_point(
      aes(x = brood_year, y = release_total, color = hatchery_location_name_fac)
    ) +
    geom_vline(aes(xintercept = cutoff_yr), linetype = "dashed") +
    facet_wrap(~treatment_control) +
    labs(
      x = "Brood year",
      y = "CWT Releases",
      color = "Hatchery"
    ) +
    scale_x_continuous(breaks = seq(2010, max_year, by = 4)) +
    scale_color_manual(values = hatch_col) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))

  plt_tog <- plt2 +
    plt3 +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "bottom")
  print(plt_tog)
  return(plt_tog)
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

soos_recs <- plot_mod_dat(mod_dat_soos, 2020, "Soos creek")
naselle_recs <- plot_mod_dat(mod_dat_naselle, 2020, "Naselle")
minter_recs <- plot_mod_dat(mod_dat_minter, 2019, "Minter creek")

# ggsave(
#   "figures/soos_recs.png",
#   soos_recs,
#   dpi = 600,
#   width = 8,
#   height = 4,
#   units = "in"
# )
# ggsave(
#   "figures/naselle_recs.png",
#   naselle_recs,
#   dpi = 600,
#   width = 8,
#   height = 4,
#   units = "in"
# )
# ggsave(
#   "figures/minter_recs.png",
#   minter_recs,
#   dpi = 600,
#   width = 8,
#   height = 4,
#   units = "in"
# )

# Model ------------------------------------------------------------------

## Soos  ------------------------------------------------------------------

#making "before" the reference category for period
mod_dat_soos$period_f <- relevel(factor(mod_dat_soos$period), ref = "Before")

#making Washington the reference category for region (has the most recoveries)
mod_dat_soos$fishery_region_f <- relevel(
  factor(mod_dat_soos$fishery_region),
  ref = "Washington"
)

#model fit
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


#diagnostics

fit2_simres <- simulateResiduals(fit2)


plot(fit2_simres)
#some (potentially minor) issues with the residuals currently

plotResiduals(fit2_simres, form = model.frame(fit2)$treatment_control)
plotResiduals(fit2_simres, form = model.frame(fit2)$fishery_region_f)
plotResiduals(fit2_simres, form = model.frame(fit2)$period_f)
testDispersion(fit2_simres)


#plotting results
plot_mod_results <- function(fit, dat) {
  EMM_2 <- emmeans(
    fit,
    ~ treatment_control * period_f * fishery_region_f,
    offset = log(1000),
    type = "response"
  )
  print(plot(EMM_2))

  #getting average offset for each group
  release_dat <- dat %>%
    select(
      brood_year,
      release_total,
      hatchery_location_name_fac,
      treatment_control,
      period
    ) %>%
    distinct()

  release_avg <- release_dat %>%
    group_by(treatment_control, period) %>%
    summarise(mean_rel = mean(release_total))

  avg_offset_rel <- dat %>%
    group_by(treatment_control, period_f, fishery_region_f) %>%
    summarise(n = n()) %>%
    left_join(release_avg, by = c("treatment_control", "period_f" = "period"))

  avg_offset_rel$period_f <- relevel(
    factor(avg_offset_rel$period_f),
    ref = "Before"
  )

  avg_offset_rel <- avg_offset_rel %>% arrange(fishery_region_f, period_f)
  print(avg_offset_rel)

  EMM_ip_2 <- emmip(
    fit,
    treatment_control ~ period_f | fishery_region_f,
    offset = log(avg_offset_rel$mean_rel),
    CIs = TRUE,
    plotit = FALSE
  )
  print(EMM_ip_2)
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

soos_plot <- plot_mod_results(fit2, mod_dat_soos) + ggtitle("a) Soos creek")
soos_plot
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


#diagnostics
fit4_simres <- simulateResiduals(fit4)


plot(fit4_simres)


testQuantiles(fit4_simres)

plotResiduals(fit4_simres, form = model.frame(fit4)$treatment_control)
plotResiduals(fit4_simres, form = model.frame(fit4)$fishery_region_f)
plotResiduals(fit4_simres, form = model.frame(fit4)$period_f)
testDispersion(fit4_simres)


#plotting results
#prints two different ways of looking at the plots
naselle_plot <- plot_mod_results(fit4, mod_dat_naselle) + ggtitle("b) Naselle")

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


#diagnostics
fit6_simres <- simulateResiduals(fit6)


plot(fit6_simres)
#some (potentially minor) issues with the residuals currently

testQuantiles(fit6_simres)

plotResiduals(fit6_simres, form = model.frame(fit6)$treatment_control)
plotResiduals(fit6_simres, form = model.frame(fit6)$fishery_region_f)

plotResiduals(fit6_simres, form = model.frame(fit6)$period_f)


#plotting results
minter_plot <- plot_mod_results(fit6, mod_dat_minter) +
  ggtitle("c) Minter creek")


#plotting all three together
plot_tog <- soos_plot /
  naselle_plot /
  minter_plot +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

plot_tog

# ggsave(
#   "figures/mod_results.png",
#   plot_tog,
#   dpi = 300,
#   width = 6,
#   height = 10,
#   units = "in"
# )

# Table for supplement ---------------------------------------------------
options(modelsummary_get = "broom")
mods <- list("Soos creek" = fit2, "Naselle" = fit4, "Minter creek" = fit6)
modelsummary(
  mods,
  shape = term ~ model + statistic,
  statistic = "conf.int",
  coef_rename = c(
    "(Intercept)" = "Intercept",
    "treatment_controlTreatment" = "Treatment",
    "period_fAfter" = "After",
    "fishery_region_fOregon" = "Oregon",
    "fishery_region_fBC" = "British Columbia",
    "fishery_region_fAlaska" = "Alaska",
    "fishery_region_fWashington" = "Washington"
  ),
  stars = TRUE,
  output = "modelsummary.docx"
)


# Random effects plot ----------------------------------------------------

re_Soos <- ranef(fit2, condVar = TRUE)
re_Soos_df <- as.data.frame(re_Soos) %>%
  add_column(trmt = "Soos creek") %>%
  mutate(brood_year = as.numeric(as.character(grp)))

re_Naselle <- ranef(fit4, condVar = TRUE)
re_Naselle_df <- as.data.frame(re_Naselle) %>%
  add_column(trmt = "Naselle") %>%
  mutate(brood_year = as.numeric(as.character(grp)))
re_Minter <- ranef(fit6, condVar = TRUE)
re_Minter_df <- as.data.frame(re_Minter) %>%
  add_column(trmt = "Minter creek") %>%
  mutate(brood_year = as.numeric(as.character(grp)))

#combine data frames for plotting
re_all <- bind_rows(re_Soos_df, re_Naselle_df, re_Minter_df)

#plots conditional mode and 2*conditional standard deviation
plt_all_re <- re_all %>%
  ggplot(aes(y = condval, x = brood_year)) +
  geom_point(aes(y = condval), size = 2) +
  geom_linerange(aes(
    ymin = condval - 2 * condsd,
    ymax = condval + 2 * condsd
  )) +
  facet_wrap(~trmt, nrow = 3) +
  xlab("Brood year \n") +
  ylab("\nConditional Mode") +
  scale_x_continuous(breaks = seq(2010, 2022, by = 2)) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white")
  )

plt_all_re

# ggsave(
#   "figures/random_effects.png",
#   plt_all_re,
#   width = 6,
#   height = 6,
#   units = "in"
# )

# Data summary for paper -------------------------------------------------
mod_dat_soos %>%
  group_by(treatment_control, period) %>%
  summarise(
    n_rec_tot = sum(n_recovered),
    avg_rec = mean(n_recovered),
    tot_rel = sum(release_total),
    avg_rel = mean(release_total)
  )
mod_dat_naselle %>%
  group_by(treatment_control, period) %>%
  summarise(
    n_rec_tot = sum(n_recovered),
    avg_rec = mean(n_recovered),
    tot_rel = sum(release_total),
    avg_rel = mean(release_total)
  )
mod_dat_minter %>%
  group_by(treatment_control, period) %>%
  summarise(
    n_rec_tot = sum(n_recovered),
    avg_rec = mean(n_recovered),
    tot_rel = sum(release_total),
    avg_rel = mean(release_total)
  )

mod_dat_soos %>%
  group_by(hatchery_location_name_fac) %>%
  summarise(
    max = max(n_recovered),
    by = brood_year[which(n_recovered == max(n_recovered))]
  )
mod_dat_naselle %>%
  group_by(hatchery_location_name_fac) %>%
  summarise(
    max = max(n_recovered),
    by = brood_year[which(n_recovered == max(n_recovered))]
  )
mod_dat_minter %>%
  group_by(hatchery_location_name_fac) %>%
  summarise(
    max = max(n_recovered),
    by = brood_year[which(n_recovered == max(n_recovered))]
  )


# Recovery rates for paper supplement-----------------------------------------------

new_dat_soos <- expand_grid(
  treatment_control = c("Control", "Treatment"),
  period = c("Before", "After"),
  fishery_region = c("BC", "Oregon", "Washington"),
  brood_year = NA,
  log_rel_total = 0
) %>%
  mutate(
    period_f = relevel(
      factor(period),
      ref = "Before"
    ),
    fishery_region_f = relevel(
      factor(fishery_region),
      ref = "Washington"
    )
  )

#get predictions
pred_soos <- predict(fit2, newdata = new_dat_soos, type = "response")
# add it to data frame
soos_tab <- new_dat_soos %>% add_column(predicted_rate = pred_soos)


new_dat_naselle <- expand_grid(
  treatment_control = c("Control", "Treatment"),
  period = c("Before", "After"),
  fishery_region = c("BC", "Alaska", "Washington"),
  brood_year = NA,
  log_rel_total = 0
) %>%
  mutate(
    period_f = relevel(
      factor(period),
      ref = "Before"
    ),
    fishery_region_f = relevel(
      factor(fishery_region),
      ref = "BC"
    )
  )

#get predictions
pred_naselle <- predict(fit4, newdata = new_dat_naselle, type = "response")
# add it to data frame
naselle_tab <- new_dat_naselle %>% add_column(predicted_rate = pred_naselle)

new_dat_minter <- expand_grid(
  treatment_control = c("Control", "Treatment"),
  period = c("Before", "After"),
  fishery_region = c("BC", "Alaska", "Washington"),
  brood_year = NA,
  log_rel_total = 0
) %>%
  mutate(
    period_f = relevel(
      factor(period),
      ref = "Before"
    ),
    fishery_region_f = relevel(
      factor(fishery_region),
      ref = "Washington"
    )
  )

#get predictions
pred_minter <- predict(fit6, newdata = new_dat_minter, type = "response")
# add it to data frame
minter_tab <- new_dat_minter %>% add_column(predicted_rate = pred_minter)

soos_tab <- soos_tab %>% add_column("Model" = "Soos creek")
naselle_tab <- naselle_tab %>% add_column("Model" = "Naselle")
minter_tab <- minter_tab %>% add_column("Model" = "Minter creek")

tabs <- bind_rows(soos_tab, naselle_tab, minter_tab) %>%
  select(-c(log_rel_total, period_f, fishery_region_f, brood_year)) %>%
  pivot_wider(
    id_cols = c(treatment_control, period, fishery_region),
    names_from = Model,
    values_from = predicted_rate
  )

tabs
#write_csv(tabs, "model_predicted_rate.csv")
