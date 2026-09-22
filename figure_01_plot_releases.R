library(dplyr)
library(ggplot2)
library(viridis)

# these are all recoveries from washington state after 1995
rel <- read.csv("data/CSV12935.TXT") |>
  dplyr::select(-record_code) |>
  dplyr::filter(brood_year >= 2010, brood_year <= 2022) |>
  dplyr::rename(tag_code = tag_code_or_release_id)
rel$cwt_1st_mark_count[which(is.na(rel$cwt_1st_mark_count))] <- 0
rel$cwt_2nd_mark_count[which(is.na(rel$cwt_2nd_mark_count))] <- 0
rel$cwt_total <- rel$cwt_1st_mark_count + rel$cwt_2nd_mark_count

rel$non_cwt_1st_mark_count[which(is.na(rel$non_cwt_1st_mark_count))] <- 0
rel$non_cwt_2nd_mark_count[which(is.na(rel$non_cwt_2nd_mark_count))] <- 0
rel$total_releases <- rel$cwt_total + rel$non_cwt_1st_mark_count + rel$non_cwt_2nd_mark_count

rel$hatchery_location_name[which(rel$hatchery_location_name == "LUMMI HATCHERY -POND")] <- "LUMMI SEA PONDS"
rel$hatchery_location_name[which(rel$hatchery_location_name == "KLICKITAT HATCHERY (YKFP)")] <- "KLICKITAT HATCHERY"

hatcheries <- data.frame(hatchery_location_name = c(
  "MARBLEMOUNT HATCHERY",
  "WELLS HATCHERY",
  "SPRING CR NFH",
  "SOOS CREEK HATCHERY",
  "BONNEVILLE HATCHERY",
  "SKOOKUM CR HATCHERY",
  "WALLACE R HATCHERY",
  "WHITE RIVER HATCHERY",
  "QUINAULT LK HATCHERY",
  "KENDALL CR HATCHERY",
  "SAMISH HATCHERY",
  "HUPP SPRINGS REARING",
  "LEWIS RIVER HATCHERY",
  "MINTER CR HATCHERY",
  "FORKS CREEK HATCHERY",
  "CLARKS CRK HATCHERY",
  "SOLDUC HATCHERY",
  "NASELLE HATCHERY",
  "BEAR SPRINGS 1  (20)",
  "WHATCOM CR HATCHERY",
  "LUMMI SEA PONDS",
  "BERNIE GOBIN HATCH",
  "WILLARD NFH",
  "LTL WHITE SALMON NFH",
  "KLICKITAT HATCHERY"
))

hatcheries$pretty_name <- c(
  "Marblemount",
  "Wells",
  "Spring Creek",
  "Soos Creek",
  "Bonneville",
  "Skookum Creek",
  "Wallace River",
  "White River",
  "Quinalt Lake",
  "Kendall Creek",
  "Samish",
  "Hupp Springs",
  "Lewis River",
  "Minter Creek",
  "Forks Creek",
  "Clarks Creek",
  "Solduc",
  "Naselle",
  "Bear Springs",
  "Whatcom Creek",
  "Lummi",
  "Bernie Gobin",
  "Willard",
  "Little White Salmon",
  "Klickitat"
)

rel <- dplyr::left_join(rel, hatcheries)

rel$release_year <- as.numeric(substr(rel$first_release_date, 1, 4))
rel$is_treatment <- ifelse(!is.na(rel$pretty_name), 1, 0)

rel_summary <- rel |>
  dplyr::group_by(release_year) |>
  dplyr::summarise(
    n_trt = sum(cwt_total[which(is_treatment == 1)]),
    n_tot = sum(cwt_total),
    n_ctrl = n_tot - n_trt,
    p_trt = n_trt / n_tot
  )

# filter out only hatcheries in the program
subset <- dplyr::filter(rel, !is.na(pretty_name))

# subset <- dplyr::filter(subset, hatchery_location_name == "SOOS CREEK HATCHERY", run==3)
# dplyr::group_by(subset, brood_year) |>
#   dplyr::summarise(annual_total = sum(total_releases)) |>
#   as.data.frame()

subset$run[which(subset$run %in% c(3, 8))] <- "Fall"
subset$run[which(subset$run %in% c(1))] <- "Spring"
subset$run[which(subset$run %in% c(2))] <- "Summer"
subset$run <- factor(subset$run, levels = c("Spring", "Summer", "Fall"))
dplyr::filter(subset, pretty_name != "Whatcom Creek") |>
  dplyr::group_by(brood_year, pretty_name, run) |>
  dplyr::summarise(n_cwt = sum(cwt_total)) |>
  ggplot(aes(brood_year, n_cwt, color = run)) +
  geom_line() +
  geom_point() +
  facet_wrap(~pretty_name, scale = "free_y", ncol = 4) +
  theme_bw() +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, big.mark = "")) +
  ylab("Coded wire tag releases (thousands)") +
  xlab("Brood year") +
  scale_x_continuous(breaks = c(2010, 2015, 2020)) +
  theme(
    # Makes the facet label background white
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 7)
  ) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave("figures/S1_release_trends_cwt.png", width = 7, height = 6)


dplyr::filter(subset, pretty_name != "Whatcom Creek", release_year <= 2023) |>
  dplyr::group_by(brood_year, pretty_name, run) |>
  dplyr::summarise(n_tot = sum(total_releases)) |>
  ggplot(aes(brood_year, n_tot, color = run)) +
  geom_line() +
  geom_point() +
  facet_wrap(~pretty_name, scale = "free_y", ncol = 4) +
  theme_bw() +
  scale_y_continuous(labels = scales::label_number(scale = 1e-3, big.mark = "")) +
  ylab("Coded wire tag releases (thousands)") +
  xlab("Brood year") +
  scale_x_continuous(breaks = c(2010, 2015, 2020)) +
  theme(
    # Makes the facet label background white
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 7)
  ) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave("figures/S2_release_trends_all.png", width = 7, height = 6)

dplyr::filter(subset, pretty_name != "Whatcom Creek", release_year <= 2023) |>
  dplyr::group_by(release_year, pretty_name, run) |>
  dplyr::summarise(
    n_cwt = sum(cwt_total),
    n_tot = sum(total_releases),
    p_cwt = n_cwt / n_tot
  ) |>
  ggplot(aes(release_year, p_cwt, color = run)) +
  geom_line() +
  geom_point() +
  facet_wrap(~pretty_name, scale = "free_y", ncol = 4) +
  theme_bw() +
  ylab("Percent of relases that are CWT") +
  xlab("Release year") +
  theme(
    # Makes the facet label background white
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 7)
  ) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave("figures/S3_release_trends_all.png", width = 7, height = 6)

# Make figure 1, showing controls and treatments

subset <- dplyr::filter(rel, hatchery_location_name %in%
  c(
    "NASELLE HATCHERY",
    "SOOS CREEK HATCHERY",
    "MINTER CR HATCHERY",
    "GORST CR REARING PND",
    "GROVERS CR HATCHERY",
    "ISSAQUAH HATCHERY",
    "VOIGHTS CR HATCHERY",
    "NEMAH HATCHERY",
    "CLEAR CREEK HATCHERY",
    "KALAMA CR HATCHERY",
    "TUMWATER FALLS HATCHERY"
  ), run == 3)

subset$release_year <- as.numeric(substr(subset$first_release_date, 1, 4))

subset$group <- NA
subset$group[which(subset$hatchery_location_name == "NASELLE HATCHERY")] <- "Naselle"
subset$group[which(subset$hatchery_location_name == "NEMAH HATCHERY")] <- "Naselle (Control)"

subset$group[which(subset$hatchery_location_name == "SOOS CREEK HATCHERY")] <- "Soos Creek"
subset$group[which(subset$hatchery_location_name %in% c(
  "GORST CR REARING PND",
  "GROVERS CR HATCHERY",
  "ISSAQUAH HATCHERY",
  "VOIGHTS CR HATCHERY"
))] <- "Soos Creek (Control)"

subset$group[which(subset$hatchery_location_name == "MINTER CR HATCHERY")] <- "Minter Creek"
subset$group[which(subset$hatchery_location_name %in% c(
  "CLEAR CREEK HATCHERY",
  "KALAMA CR HATCHERY",
  "TUMWATER FALLS HATCHERY"
))] <- "Minter Creek (Control)"

subset |>
  dplyr::group_by(group, release_year) |>
  dplyr::summarise(n_tot = sum(total_releases)) |>
  ggplot(aes(release_year, n_tot)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~group, scale = "free_y", ncol = 2) +
  theme_bw() +
  ylab("All releases") +
  xlab("Release year") +
  theme(
    # Makes the facet label background white
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 7)
  ) +
  scale_color_viridis()
ggsave("figures/Figure1_release_trends_control_treatment.png", width = 7, height = 6)
