# species with data from wildlife portal
portal_spp <- read.csv("E:/NFWF_PhaseII/Alaska/Data/AK_wildlife_portal/wildlife_portal_species.csv")

# species currently in the assessment
assessment_spp <- read.csv("E:/NFWF_PhaseII/Alaska/Data/AK_spp_list_11-18-21.csv")
assessment_spp$include <- 'in_assessment'

# species in the SWAP
swap_spp <- read.csv("E:/NFWF_PhaseII/Alaska/Data/AK_SWAP_species.csv")
swap_spp$SWAP <- 'in SWAP'

spp_merge <- merge(
  portal_spp, swap_spp, by.x='scientific_name', by.y='Scientific.Name',
  all=TRUE)
in_swap_not_portal <- spp_merge[is.na(spp_merge$elcode), ]
in_portal_not_swap <- spp_merge[is.na(spp_merge$SWAP), ]

in_swap_and_portal <- merge(
  portal_spp, swap_spp, by.x='scientific_name', by.y='Scientific.Name')
add_test <- merge(
  in_swap_and_portal, assessment_spp, by.x='scientific_name', by.y='scientific',
  all=TRUE)
add_spp <- add_test[is.na(add_test$include), ]

# subsistence species list
non_fish_subs_resources <- read.csv(
  "E:/Datasets/AK_Data_Portal/Subsistence_Harvests%3A_Non-Fishing_Resources-shp/Subsistence_harvests_nonfishing_resources_count_by_spp.csv")
fish_subs_resources <- read.csv(
  "E:/Datasets/AK_Data_Portal/Subsistence_Harvests%3A_Fishing_Resources-shp/Subsistence_harvest_fishing_resources_count_by_spp.csv")
subsistence_spp <- rbind(non_fish_subs_resources, fish_subs_resources)
subsistence_spp <- subsistence_spp[, c('Resource', 'FREQUENCY')]
write.csv(subsistence_spp, "E:/NFWF_PhaseII/Alaska/Data/Subsistence_spp/AK_Data_Portal_subsistence_harvests_spp_count_combined.csv")

# subsistence spp with notes
subsistence_spp <- read.csv("E:/NFWF_PhaseII/Alaska/Data/Subsistence_spp/AK_Data_Portal_subsistence_harvests_spp_count_combined.csv")
