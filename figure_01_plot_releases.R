library(dplyr)

# these are all recoveries from washington state after 1995
rel <- read.csv("data/releases_2010_2026.csv") |>
  dplyr::select(-record_code) |>
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

# filter out only hatcheries in the program
subset <- dplyr::filter(rel, !is.na(pretty_name))

subset$release_year <- as.numeric(substr(subset$first_release_date, 1, 4))

subset$run[which(subset$run %in% c(3, 8))] <- "Fall"
subset$run[which(subset$run %in% c(1))] <- "Spring"
subset$run[which(subset$run %in% c(2))] <- "Summer"
subset$run <- factor(subset$run, levels = c("Spring", "Summer", "Fall"))
dplyr::filter(subset, pretty_name != "Whatcom Creek", release_year <= 2023) |>
  dplyr::group_by(release_year, pretty_name, run) |>
  dplyr::summarise(n_cwt = sum(cwt_total)) |>
  ggplot(aes(release_year, n_cwt, color = run)) +
  geom_line() +
  geom_point() +
  facet_wrap(~pretty_name, scale = "free_y", ncol = 4) +
  theme_bw() +
  ylab("Coded wire tag releases") +
  xlab("Release year") +
  theme(
    # Makes the facet label background white
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 7)
  ) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave("figures/01_release_trends.png", width = 7, height = 6)


dplyr::filter(subset, pretty_name != "Whatcom Creek", release_year <= 2023) |>
  dplyr::group_by(release_year, pretty_name, run) |>
  dplyr::summarise(n_tot = sum(total_releases)) |>
  ggplot(aes(release_year, n_tot, color = run)) +
  geom_line() +
  geom_point() +
  facet_wrap(~pretty_name, scale = "free_y", ncol = 4) +
  theme_bw() +
  ylab("All releases") +
  xlab("Release year") +
  theme(
    # Makes the facet label background white
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 7)
  ) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8)
ggsave("figures/S1_release_trends_all.png", width = 7, height = 6)
