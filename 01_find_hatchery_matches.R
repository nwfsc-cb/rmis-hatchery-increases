library(dplyr)
library(tidyr)

# these are all recoveries from washington state after 1995
d <- readRDS("wa_data_after_1995.rds")

# let's only look at ocean commercial fisheries
d <- dplyr::filter(
  d, fishery == 10,
  !is.na(recovery_rmis_region),
  recovery_rmis_region != "",
  brood_year >= 2010
)

# group by hatchery and brood year and summarize effective number
n_hatch_by <- dplyr::group_by(d, hatchery_location_name, brood_year, recovery_rmis_region) |>
  dplyr::summarize(
    est_num = sum(estimated_number, na.rm = T),
    .groups = "drop_last"
  )


# Calculate proportions within each brood year for each hatchery
hatch_props_yearly <- n_hatch_by |>
  dplyr::group_by(hatchery_location_name, brood_year) |>
  dplyr::mutate(prop = est_num / sum(est_num, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(prop = ifelse(is.nan(prop), 0, prop))

# Create the wide matrix, keeping brood_year as a row identifier
matrix_data_yearly <- hatch_props_yearly |>
  dplyr::select(hatchery_location_name, brood_year, recovery_rmis_region, prop) |>
  tidyr::pivot_wider(
    names_from = recovery_rmis_region,
    values_from = prop,
    values_fill = 0
  )

# Isolate Soos Creek's yearly distribution
soos_yearly <- matrix_data_yearly |>
  dplyr::filter(hatchery_location_name == this_hatchery)

# Join other hatcheries to Soos Creek on the same brood year
# and calculate the Hellinger distance for each year
# Hellinger is sqrt transforming then (1/2) * sum((sqrt(other_vec) - sqrt(soos_vec))^2)
yearly_distances <- matrix_data_yearly |>
  dplyr::filter(hatchery_location_name != this_hatchery) |>
  dplyr::inner_join(soos_yearly, by = "brood_year", suffix = c("_other", "_soos")) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    h_distance = {
      # where(is.numeric) forces it to ignore the character hatchery name columns!
      other_vec <- c_across(c(ends_with("_other") & where(is.numeric)))
      soos_vec <- c_across(c(ends_with("_soos") & where(is.numeric)))

      # Hellinger formula
      sqrt(0.5 * sum((sqrt(other_vec) - sqrt(soos_vec))^2))
    }
  ) |>
  dplyr::ungroup()

# Find the overall best match across all shared years
# (We average the distances across all years they have in common)
final_ranking <- yearly_distances |>
  dplyr::group_by(hatchery_location_name_other) |>
  dplyr::summarize(
    mean_hellinger_distance = mean(h_distance),
    years_compared = n()
  ) |>
  dplyr::arrange(mean_hellinger_distance) |>
  dplyr::rename(hatchery_location_name = hatchery_location_name_other)

final_ranking <- dplyr::filter(final_ranking, years_compared >= 10)

#> final_ranking
# A tibble: 39 × 3
# hatchery_location_name mean_hellinger_distance years_compared
# <chr>                                    <dbl>          <int>
#   1 GEORGE ADAMS HATCHERY                    0.324             13
# 2 GROVERS CR HATCHERY                      0.334             13
# 3 HOODSPORT HATCHERY                       0.336             13
# 4 CLEAR CREEK HATCHERY                     0.372             13
# 5 MINTER CR HATCHERY                       0.384             10
# 6 HUPP SPRINGS REARING                     0.404             10
# 7 SAMISH HATCHERY                          0.420             12
# 8 SKOOKUM CR HATCHERY                      0.434             11
# 9 WHITE RIVER HATCHERY                     0.453             10
# 10 LYONS FERRY HATCHERY                     0.456             13
#
