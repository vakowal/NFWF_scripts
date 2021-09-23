# location where NWI categories should be saved for each IUCN habitat class
intermediate_dir <- "D:/NFWF_PhaseIII/Great_Lakes/Data - wildlife/NWI/IUCN_match_tables"

# select NWI rows
subset_list <- list()
nwi_df <- read.csv("D:/Datasets/FWS/NWI/NWI-Code-Definitions-Oct-19/NWI_Code_Definitions/NWI_Code_Definitions_GK_rearrange.csv")

##################################################################
# Estuarine systems (tidal wetlands with sporadic access to ocean)
# skip for Great Lakes ???
estuar_sub <- nwi_df[nwi_df$SYSTEM_NAM == 'Estuarine', ]

# 13.4, coastal brackish/saline lagoons or marine lakes
# 12.7, mangrove submerged roots
# 9.1, estuaries

##################################
# Lacustrine systems (lacus, lake)
lacus_sub <- nwi_df[nwi_df$SYSTEM_NAM == 'Lacustrine', ]

# 15.3, aquaculture ponds
# 15.2, artificial ponds
# 15.1, artificial water storage areas
aqua <- lacus_sub[((lacus_sub$FIRST_MODI == 'Managed') | (lacus_sub$FIRST_MODI == 'Excavated')) &
                    (lacus_sub$WATER_REGI == 'Permanently Flooded'), ]
aqua$IUCN_class <- '15.3'
subset_list[['15.3']] <- aqua
aqua1 <- aqua
aqua1$IUCN_class <- '15.2'
subset_list[['15.2']] <- aqua1
aqua2 <- aqua
aqua2$IUCN_class <- '15.1'
subset_list[['15.1']] <- aqua2

# 13.5, coastal freshwater lakes  # skip for Great Lakes

# 5.5, permanent freshwater lakes
# 5.7, permanent freshwater marshes/pools under 8ha
subs1 <- lacus_sub[(lacus_sub$WATER_REGI == 'Permanently Flooded') & (lacus_sub$FIRST_MODI == ""), ]
subs1$IUCN_class <- '5.5'
subset_list[['5.5']] <- subs1
subs2 <- subs1
subs2$IUCN_class <- '5.7'
subset_list[['5.7']] <- subs2

# 5.6, seasonal/intermittent freshwater lakes over 8ha
# 5.8, seasonal/intermittent freshwater marshes/pools under 8ha
subs1 <- lacus_sub[((lacus_sub$WATER_REGI == 'Seasonally Flooded') | (lacus_sub$WATER_REGI == 'Intermittently Flooded')) &
                     (lacus_sub$FIRST_MODI == ""), ]
subs1$IUCN_class <- '5.6'
subset_list[['5.6']] <- subs1
subs2 <- subs1
subs2$IUCN_class <- '5.8'
subset_list[['5.8']] <- subs2

# 5.14, saline, brackish, alkaline lakes
saline_subs <- subset(lacus_sub, grepl("saline", lacus_sub$FIRST_MODI))
saline_subs1 <- saline_subs[saline_subs$WATER_RE_1 == 'Nontidal', ]
alkaline_subs <- subset(lacus_sub, grepl("Alkaline", lacus_sub$FIRST_MODI))
union <- rbind(saline_subs1, alkaline_subs)
union$IUCN_class <- '5.14'
subset_list[['5.14']] <- union

################
# Marine systems : skip these for Great Lakes ???
marine_sub <- nwi_df[nwi_df$SYSTEM_NAM == 'Marine', ]

# 11.1, marine deep benthic continental slope
# 12, marine intertidal
# 12.4, marine intertidal: mud flats and salt flats
# 12.5, marine intertidal: salt marshes
# 9.8, marine neritic: coral reef


###################################
# Palustrine systems (palus, marsh)
palus_sub <- nwi_df[nwi_df$SYSTEM_NAM == 'Palustrine', ]

# 5.16, permanent saline, brackish, alkaline inland marshes
saline_subs <- subset(palus_sub, grepl("saline", palus_sub$FIRST_MODI))
saline_subs1 <- saline_subs[saline_subs$WATER_RE_1 == 'Nontidal', ]
alkaline_subs <- subset(palus_sub, grepl("Alkaline", palus_sub$FIRST_MODI))
union <- rbind(saline_subs1, alkaline_subs)
union$IUCN_class <- '5.16'
subset_list[['5.16']] <- union

# 5.3, shrub dominated wetlands
shrub_1 <- subset(palus_sub, grepl("Scrub-Shrub", palus_sub$CLASS_NAME))
shrub_2 <- shrub_1[(shrub_1$WATER_RE_1 == 'Nontidal'), ]
shrub_3 <- shrub_2[shrub_2$FIRST_MODI == "", ]
shrub_3$IUCN_class <- '5.3'
subset_list[['5.3']] <- shrub_3

# 5.4, bogs, marshes, swamps, fens, peatlands
attributes_5.4 <- c('PML/EM1B', 'PML1/2B', 'PABFg', 'PEM1/AB3Fg', 'PEM1/ABFg', 'PEM1/ML1B', 'PEM1/ML1C', 'PEM1/ML1E', 'PEM1/ML2B', 'PEM1/ML2C', 'PEM1/MLB', 'PML/EM1B', 'PML1/2B', 'PML1/2C', 'PML1/EM1B', 'PML1/EM1C', 'PML1B', 'PML1C', 'PML1E', 'PML2/1B', 'PML2/EM1B', 'PML2B', 'PMLB')
bogs <- subset(palus_sub, palus_sub$ATTRIBUTE %in% attributes_5.4)
bogs$IUCN_class <- '5.4'
subset_list[['5.4']] <- bogs

##################
# Riverine systems

river_sub <- nwi_df[nwi_df$SYSTEM_NAM == 'Riverine', ]

# 15.9, artificial canals and drainage ditches
subs1 <- river_sub[((river_sub$FIRST_MODI == 'Managed') | (river_sub$FIRST_MODI == 'Excavated')) &
                     (river_sub$WATER_REGI == 'Permanently Flooded'), ]
subs1$IUCN_class <- '15.9'
subset_list[['15.9']] <- subs1

# 5.9, inland freshwater springs and oases
# no matches in NWI

# 5.1, permanent rivers/streams/creeks
# 5.13, permanent inland deltas
subs1 <- river_sub[(river_sub$FIRST_MODI == "") &
                     (river_sub$WATER_REGI == 'Permanently Flooded'), ]
subs1$IUCN_class <- '5.1'
subset_list[['5.1']] <- subs1
subs2 <- subs1
subs2$IUCN_class <- '5.13'
subset_list[['5.13']] <- subs2

# 5.18, karst and other subterranean hydrological systems
# no matches in NWI

# 5.2, seasonal/intermittent rivers, streams, creeks
subs1 <- river_sub[(river_sub$FIRST_MODI == "") &
                     ((river_sub$WATER_REGI == 'Seasonally Flooded') |
                      (river_sub$WATER_REGI == 'Intermittently Flooded')), ]
subs1$IUCN_class <- '5.2'
subset_list[['5.2']] <- subs1

# rbind all subsets together
matched_NWI <- do.call(rbind, subset_list)
matched_table <- matched_NWI[, c('ATTRIBUTE', 'IUCN_class')]
write.csv(
  matched_table, paste(intermediate_dir, 'IUCN_NWI_match.csv', sep='/'),
  row.names=FALSE)
