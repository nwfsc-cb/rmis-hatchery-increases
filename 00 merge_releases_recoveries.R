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

# these are all recoveries from washington state after 1995
rel <- read.csv("data/releases_2010_2026.csv") |>
  dplyr::select(-record_code) |>
  dplyr::rename(tag_code = tag_code_or_release_id)
rel$cwt_1st_mark_count[which(is.na(rel$cwt_1st_mark_count))] <- 0
rel$cwt_2nd_mark_count[which(is.na(rel$cwt_2nd_mark_count))] <- 0
rel$cwt_total <- rel$cwt_1st_mark_count + rel$cwt_2nd_mark_count

# tags_to_keep <- which(recoveries$tag_code %in% rel$tag_code)
# recoveries <- recoveries[tags_to_keep,]
# filter out fish that were released outside of OR or WA -- these
# releases are just from those states. this also includes potentially
# hatcheries that released fish but
joined <- dplyr::left_join(recoveries, rel, by = "tag_code") |>
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
  recovery_location_name,
  stock_location_name,
  release_location_state,
  release_location_rmis_region,
  release_location_rmis_basin,
  recovery_location_code
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

write.csv(joined, "data/joined_data.csv")
