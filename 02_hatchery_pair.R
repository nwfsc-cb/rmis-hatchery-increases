# specify the treatment hatchery
this_hatchery <- "SOOS CREEK HATCHERY"

source("01_find_hatchery_matches.R")
control_hatchery <- final_ranking$hatchery_location_name[1]

# now filter the recoveries to be just from these hatcheries
# summarize by brood year and recovery region
# though we could do this at an even coarser scale - BACI would be before / after 2021?
recoveries <- readRDS("wa_data_after_1995.rds") |>
  dplyr::filter(hatchery_location_name %in% c(this_hatchery, control_hatchery)) |>
  dplyr::filter(brood_year > 2000) |>
  dplyr::group_by(hatchery_location_name, brood_year, recovery_rmis_region) |>
  dplyr::summarise(
    n = n(),
    n_est = sum(estimated_number, na.rm = TRUE),
    .groups = "drop_last"
  )

releases <- readRDS("all_releases_june2026.rds")
releases$cwt_1st_mark_count[which(is.na(releases$cwt_1st_mark_count))] <- 0
releases$cwt_2nd_mark_count[which(is.na(releases$cwt_2nd_mark_count))] <- 0
releases$cwt <- releases$cwt_1st_mark_count + releases$cwt_2nd_mark_count
releases <- dplyr::filter(releases, hatchery_location_name %in% c(this_hatchery, control_hatchery)) |>
  dplyr::filter(brood_year >= min(recoveries$brood_year)) |>
  dplyr::group_by(brood_year, hatchery_location_name) |>
  dplyr::summarise(n_cwt_release = sum(cwt, na.rm = T), .groups = "drop_last")
