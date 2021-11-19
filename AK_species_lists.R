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
