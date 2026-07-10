joined <- read.csv("data/joined_data.csv")

data_list <- list()

# This is probably easiest to do for each treatment hatchery


# These are fall
sub <- dplyr::filter(joined, run == 3, hatchery_location_name %in% c("CLARKS CRK HATCHERY", "SOOS CREEK HATCHERY", "GROVERS CR HATCHERY", "GORST CR REARING PND", "ISSAQUAH HATCHERY", "VOIGHTS CR HATCHERY"))
sub$category <- ifelse(sub$hatchery_location_name %in% c("CLARKS CRK HATCHERY", "SOOS CREEK HATCHERY"), "treatment", "control")
data_list[[1]] <- sub

# These are fall
sub <- dplyr::filter(joined, run == 3, hatchery_location_name %in% c("QUINAULT LK HATCHERY", "SALMON R FISH CULTUR"))
sub$category <- ifelse(sub$hatchery_location_name == "QUINAULT LK HATCHERY", "treatment", "control")
data_list[[2]] <- sub

# These are fall
sub <- dplyr::filter(joined, run == 3, hatchery_location_name %in% c("SAMISH HATCHERY", "GLENWOOD SPRINGS"))
sub$category <- ifelse(sub$hatchery_location_name %in% c("SAMISH HATCHERY"), "treatment", "control")
data_list[[3]] <- sub

# These are spring
sub <- dplyr::filter(joined, run == 1, hatchery_location_name %in% c("LEWIS RIVER HATCHERY", "COWLITZ SALMON HATCHERY"))
sub$category <- ifelse(sub$hatchery_location_name %in% c("LEWIS RIVER HATCHERY"), "treatment", "control")
data_list[[4]] <- sub

# These are fall
sub <- dplyr::filter(joined, run == 3, hatchery_location_name %in% c("FORKS CREEK HATCHERY", "NASELLE HATCHERY", "NEMAH HATCHERY"))
sub$category <- ifelse(sub$hatchery_location_name %in% c("FORKS CREEK HATCHERY", "NASELLE HATCHERY"), "treatment", "control")
data_list[[5]] <- sub

# These are fall
sub <- dplyr::filter(joined, run == 3, hatchery_location_name %in% c("MINTER CR HATCHERY", "CLEAR CREEK HATCHERY", "KALAMA CR HATCHERY", "TUMWATER FALLS HATCHERY"))
sub$category <- ifelse(sub$hatchery_location_name %in% c("MINTER CR HATCHERY"), "treatment", "control")
data_list[[6]] <- sub

saveRDS(data_list, "data/data_for_modeling.rds")
