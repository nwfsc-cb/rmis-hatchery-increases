library(dplyr)

# these are all recoveries from washington state after 1995
rel <- read.csv("data/releases_2010_2026.csv") |>
  dplyr::select(-record_code) |>
  dplyr::rename(tag_code = tag_code_or_release_id)
rel$cwt_1st_mark_count[which(is.na(rel$cwt_1st_mark_count))] <- 0
rel$cwt_2nd_mark_count[which(is.na(rel$cwt_2nd_mark_count))] <- 0
rel$cwt_total <- rel$cwt_1st_mark_count + rel$cwt_2nd_mark_count

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
  "WHATCOM CR HATCHERY"
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
  "Whatcom Creek"
)

rel <- dplyr::left_join(rel, hatcheries)

# filter out only hatcheries in the program
subset <- dplyr::filter(rel, !is.na(pretty_name))

subset$release_year <- as.numeric(substr(subset$first_release_date, 1, 4))

dplyr::filter(subset, pretty_name != "Whatcom Creek", release_year <= 2025) |>
  dplyr::group_by(release_year, pretty_name) |>
  dplyr::summarise(n_cwt = sum(cwt_total)) |>
  ggplot(aes(release_year, n_cwt)) +
  geom_line() +
  facet_wrap(~pretty_name, scale = "free_y") +
  theme_bw() +
  ylab("Coded wire tag releases") +
  xlab("Release year") +
  theme(
    # Makes the facet label background white
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 7)
  )
ggsave("figures/01_release_trends.png", width = 7, height = 6)
