# select NWI rows
nwi_df <- read.csv("D:/Datasets/FWS/NWI/NWI-Code-Definitions-Oct-19/NWI_Code_Definitions/NWI_Code_Definitions_GK_rearrange.csv")

saline_subs <- subset(nwi_df, grepl("saline", nwi_df$FIRST_MODI))
saline_subs1 <- saline_subs[saline_subs$WATER_RE_1 == 'Nontidal', ]
alkaline_subs <- subset(nwi_df, grepl("Alkaline", nwi_df$FIRST_MODI))
union <- paste(c(saline_subs1$ATTRIBUTE, alkaline_subs$ATTRIBUTE), sep=', ', collapse = ', ')

# shrub dominated wetlands
Scrub-Shrub
shrub_1 <- subset(nwi_df, grepl("saline", nwi_df$FIRST_MODI))