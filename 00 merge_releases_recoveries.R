library(dplyr)
library(tidyr)

# read in recoveries
for (y in 2010:2025) {
  f <- read.csv(paste0("data/recoveries_", y, ".csv"))
  if (y == 2010) {
    recoveries <- f
  } else {
    recoveries <- rbind(recoveries, f)
  }
}

# Load and clean release data
rel <- read.csv("data/CSV12935.TXT") |>
  dplyr::select(-record_code) |>
  dplyr::rename(tag_code = tag_code_or_release_id) |>
  mutate(
    cwt_1st_mark_count = ifelse(is.na(cwt_1st_mark_count), 0, cwt_1st_mark_count),
    cwt_2nd_mark_count = ifelse(is.na(cwt_2nd_mark_count), 0, cwt_2nd_mark_count),
    cwt_total = cwt_1st_mark_count + cwt_2nd_mark_count,
    non_cwt_1st_mark_count = ifelse(is.na(non_cwt_1st_mark_count), 0, non_cwt_1st_mark_count),
    non_cwt_2nd_mark_count = ifelse(is.na(non_cwt_2nd_mark_count), 0, non_cwt_2nd_mark_count),
    total_tags = cwt_total + non_cwt_1st_mark_count + non_cwt_2nd_mark_count
  ) |>
  filter(release_stage != "Y") |>
  mutate(release_year = as.numeric(substr(first_release_date, 1, 4)))

# Define Treatment vs. Control groups
rel <- rel |>
  mutate(
    treatment_control = case_when(
      # Treatment Groups
      release_location_name == "NASELLE R    24.0543" ~ "Treatment",
      release_location_name %in% c("MINTER CR    15.0048", "MINTER CR TR 15.0051") ~ "Treatment",
      release_location_rmis_basin == "DUWA" &
        stock_location_name %in% c("BIG SOOS CR  09.0072", "SOOS CR+GREEN R NATIVE") ~ "Treatment",
      # Control Groups
      release_location_name %in% c("NORTH NEMAH R24.0460", "N NEMAH R TR 24.0498", "M NEMAH R TR 24.0508") ~ "Control",
      release_location_name %in% c("CLEAR CR    11.0013C", "KALAMA CR    11.0017", "MCALLISTER SPRINGS HATCH") ~ "Control",
      release_location_rmis_basin == "EKPN" &
        release_location_name %in% c("GROVERS CR   15.0299", "GROVERS CR HATCHERY") ~ "Control",
      TRUE ~ NA_character_
    )
  )

# Define groupings to apply "Before" vs "After" universally across Treatment & Control
# - 1 (Naselle / Nemah): After cutoff = 2020
# - 2 (Minter / Clear / Kalama / McAllister): After cutoff = 2019
# - 3 (Duwamish / Grovers): After cutoff = 2020

rel <- rel |>
  filter(!is.na(treatment_control)) |>
  mutate(
    region = case_when(
      release_location_name %in% c("NASELLE R    24.0543", "NORTH NEMAH R24.0460", "N NEMAH R TR 24.0498", "M NEMAH R TR 24.0508") ~ "Willapa",
      release_location_name %in% c("MINTER CR    15.0048", "MINTER CR TR 15.0051", "CLEAR CR    11.0013C", "KALAMA CR    11.0017", "MCALLISTER SPRINGS HATCH") ~ "South Sound",
      release_location_rmis_basin %in% c("DUWA", "EKPN") ~ "Central Sound",
      TRUE ~ "Other"
    ),
    period = case_when(
      region == "Willapa" & release_year >= 2020 ~ "After",
      region == "South Sound" & release_year >= 2019 ~ "After",
      region == "Central Sound" & release_year >= 2020 ~ "After",
      TRUE ~ "Before"
    )
  )

# Check the BACI table
table(rel$treatment_control, rel$period)

# Keep ONLY the defined BACI study groups
rel <- dplyr::filter(rel, !is.na(treatment_control))

avg_releases <- rel |>
  filter(!is.na(treatment_control)) |>
  group_by(treatment_control, period, release_location_name) |>
  summarize(
    mean_releases = mean(cwt_total, na.rm = TRUE),
    median_releases = median(cwt_total, na.rm = TRUE),
    total_releases = sum(cwt_total, na.rm = TRUE),
    num_release_groups = n(),
    .groups = "drop"
  ) |>
  arrange(treatment_control, release_location_name, period)


# tags_to_keep <- which(recoveries$tag_code %in% rel$tag_code)
# recoveries <- recoveries[tags_to_keep,]
# filter out fish that were released outside of OR or WA -- these
# releases are just from those states. this also includes potentially
# hatcheries that released fish but
joined <- dplyr::left_join(dplyr::select(recoveries, -period), rel, by = "tag_code") |>
  dplyr::filter(
    !is.na(number_cwt_estimated),
    !is.na(cwt_total)
  )

joined <- dplyr::select(
  joined,
  tag_code, # RMIS tag code
  run_year,
  recovery_date,
  fishery,
  number_cwt_estimated, # this is the expanded estimate accounting for sampling
  related_group_id,
  run,
  brood_year,
  cwt_total,
  total_tags,
  recovery_location_name,
  stock_location_name,
  release_location_state,
  release_location_rmis_region,
  release_location_rmis_basin,
  recovery_location_code,
  hatchery_location_name,
  hatchery_location_code,
  first_release_date,
  release_stage,
  comments,
  avg_length,
  avg_weight,
  period,
  treatment_control
)

# locations
locs <- read.csv("data/rmis_locations.csv") |>
  dplyr::select(location_code, rmis_region) |>
  dplyr::rename(
    recovery_location_code = location_code,
    recovery_rmis_region = rmis_region
  ) |>
  dplyr::group_by(recovery_location_code) |>
  dplyr::summarise(recovery_rmis_region = recovery_rmis_region[1]) |>
  as.data.frame()

joined <- dplyr::left_join(joined, locs)
# release year is usually brood year + 1
# joined$release_year <- as.numeric(as.character(substr(joined$first_release_date, 1, 4)))

# add fishery regions
# joined$recovery_location_code
joined$fishery_region <- "Other"
joined$fishery_region[which(substr(joined$recovery_location_code, 1, 2) == "1M")] <- "Alaska"
joined$fishery_region[which(substr(joined$recovery_location_code, 1, 2) == "2M")] <- "BC"
joined$fishery_region[which(substr(joined$recovery_location_code, 1, 2) == "3M")] <- "Washington"
joined$fishery_region[which(substr(joined$recovery_location_code, 1, 2) == "5M")] <- "Oregon"
joined$fishery_region[which(substr(joined$recovery_location_code, 1, 2) == "6M")] <- "California"

# rel |>
#   dplyr::filter(hatchery_location_name == "CLARKS CRK HATCHERY") |>
#   dplyr::group_by(brood_year) |>
#   dplyr::summarise(n = sum(cwt_total))


write.csv(joined, "data/joined_data.csv")

##################################################################
# The whole point of this block is to do some labeling and for
# each of the treatment hatcheries, identify the RMIS region and basin
# they come from
##################################################################
treatment_hatcheries <- read.csv("data/treatment_hatcheries.csv")
treatment_hatcheries$release_location_rmis_region <- ""
treatment_hatcheries$release_location_rmis_basin <- ""
treatment_hatcheries$comment <- ""

for (i in 1:nrow(treatment_hatcheries)) {
  sub <- dplyr::filter(rel, hatchery_location_name == treatment_hatcheries$release_location_name[i])

  if (treatment_hatcheries$release_location_name[i] == "MARBLEMOUNT HATCHERY") {
    sub$release_location_rmis_basin <- "UPSK"
    treatment_hatcheries$comment[i] <- "All releases from SKAG-UPSK are from this hatchery"
  }
  if (treatment_hatcheries$release_location_name[i] == "WELLS HATCHERY") {
    sub <- dplyr::filter(sub, release_location_rmis_region == "CRGN")
    treatment_hatcheries$comment[i] <- "Mis releases in 2017-2018 from UPCR filtered out"
  }
  if (treatment_hatcheries$release_location_name[i] == "BONNEVILLE HATCHERY") {
    sub <- dplyr::filter(sub, release_location_name == "TANNER CR (BNVILLE)")
    treatment_hatcheries$comment[i] <- "Mis releases in Umatilla / Yakima filtered out"
  }
  if (treatment_hatcheries$release_location_name[i] == "MINTER CR HATCHERY") {
    sub <- dplyr::filter(sub, release_location_rmis_basin == "EKPS")
    treatment_hatcheries$comment[i] <- "Mis releases in DES filtered out"
  }
  if (treatment_hatcheries$release_location_name[i] == "CLARKS CRK HATCHERY") {
    sub <- dplyr::filter(sub, release_location_rmis_region == "MPS")
    treatment_hatcheries$comment[i] <- "Mis releases filtered out"
  }

  if (length(unique(sub$release_location_rmis_region)) == 1) {
    treatment_hatcheries$release_location_rmis_region[i] <- sub$release_location_rmis_region[1]
  }
  if (length(unique(sub$release_location_rmis_basin)) == 1) {
    treatment_hatcheries$release_location_rmis_basin[i] <- sub$release_location_rmis_basin[1]
  }
}
write.csv(treatment_hatcheries, "data/treatment_hatcheries_rmis.csv")


##################################################################
# Finally, we can try to assign some control hatcheries as hatcheries
# coming from the same RMIS basin / region. There's only 20 treatment
# hatcheries here
##################################################################
treatment_hatcheries <- read.csv("data/treatment_hatcheries_rmis.csv")

# let's use region (coarser) instead of basin (fine scale)
# for (i in 1:nrow(treatment_hatcheries)) {
#   sub <- dplyr::filter(
#     rel,
#     release_location_rmis_region == treatment_hatcheries$release_location_rmis_region[i]
#   )
# }

# NOTES:
# MARBLEMOUNT has no good controls
# WELLS has no good controls
# SPRING CR NFH -> All of the other hatcheries in the RMIS region release fish at much larger sizes, so no good controls
# SOOS CREEK -> "CLARKS CRK HATCHERY" /  "SOOS CREEK HATCHERY" treatment vs "GROVERS CR HATCHERY", "GORST CR REARING PND", "ISSAQUAH HATCHERY", "VOIGHTS CR HATCHERY" controls
# BONNEVILLE HATCHERY -> has no good controls
# SKOOKUM CR HATCHERY -> KENDALL / SKOOKUM are the only spring releases in RMIS NOWA region, no good controls
# WALLACE R HATCHERY -> not many summer equivalents
# WHITE RIVER HATCHERY -> "CLARKS CRK HATCHERY" /  "SOOS CREEK HATCHERY" treatment vs "GROVERS CR HATCHERY", "GORST CR REARING PND", "ISSAQUAH HATCHERY", "VOIGHTS CR HATCHERY" controls
# QUINAULT LK HATCHERY -> SALMON R FISH CULTUR
# KENDALL CR HATCHERY -> KENDALL / SKOOKUM are the only spring releases in RMIS NOWA region, no good controls
# SAMISH HATCHERY -> release stage isn't indicated, and KENDALL / SKOOKUM CR HATCHERY / SAMISH HATCHERY all treatments are are the only hatcheries in NOWA (other than Glenwood springs)
# HUPP SPRINGS REARING -> there are no spring equivalents in RMIS SPS region
# LEWIS RIVER HATCHERY -> BIG CR HATCHERY, CEDC NET PENS,CLACKAMAS HATCHERY, COWLITZ SALMON HATCHERY, DEEP R NET PENS, DEXTER PONDS (WILLAM, FALLERT CR HATCHERY, GNAT CR HATCHERY, KLASKANINE HATCHERY,LEABURG HATCHERY,MCKENZIE HATCHERY,MINTO PONDS (N SANTIAM R),SANDY HATCHERY,SOUTH SANTIAM HATCH,SPEELYAI HATCHERY,WILLAMETTE HATCHERY
# MINTER CR HATCHERY ->  "CLEAR CREEK HATCHERY", "KALAMA CR HATCHERY", "TUMWATER FALLS HATCHERY" are controls
# FORKS CREEK HATCHERY -> lumped with NASELLE, Nemah as control
# CLARKS CRK HATCHERY -> CLARKS CRK HATCHERY / SOOS CREEK HATCHERY / WHITE RIVER HATCHERY are all in RMIS REGION MPS and are treatments, with GROVERS CR HATCHERY / GORST CR REARING PND / ISSAQUAH HATCHERY / VOIGHTS CR HATCHERY being controls (at least for 2021)
# SOLDUC HATCHERY -> not enough summer equivalents in this region
# NASELLE HATCHERY -> lumped with Forks Creek, Nemah as control
# BEAR SPRINGS 1  (20) -> not enough brood years
# WHATCOM CR HATCHERY -> release stage isn't indicated, 1 release
# BERNIE GOBIN -> these are summer run, and the major hatchery in NPS producing summer run (= Wallace) is also part of the program
# WILLARD -> these are late fall (run = 8) and the only other CECR hatcheries producing these (Klickitat, Ltl White Salmon) are part of the program
