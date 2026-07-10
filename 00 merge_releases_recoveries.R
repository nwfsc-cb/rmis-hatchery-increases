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

# delete releases that have comments indicating vaccination, trucking, etc
# keeps 2948 / 3880 groups
bad_comments <- read.csv("data/bad_comments.csv")
rel <- dplyr::filter(rel, comments %in% bad_comments$comments == FALSE)

# also filter out fry / fingerling / pre-smolts rlease stages
# rel <- dplyr::filter(rel, release_stage %in% c("S", "Y"))

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
  recovery_location_code,
  hatchery_location_name,
  hatchery_location_code,
  first_release_date,
  release_stage,
  comments,
  avg_length,
  avg_weight
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
joined$release_year <- as.numeric(as.character(substr(joined$first_release_date, 1, 4)))

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
#
# Finally, we can try to assign some control hatcheries as hatcheries
# coming from the same RMIS basin / region. There's only 20 treatment
# hatcheries here
#
##################################################################
treatment_hatcheries <- read.csv("data/treatment_hatcheries_rmis.csv")

# let's use region (coarser) instead of basin (fine scale)
for (i in 1:nrow(treatment_hatcheries)) {
  sub <- dplyr::filter(
    rel,
    release_location_rmis_region == treatment_hatcheries$release_location_rmis_region[i]
  )
}

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
