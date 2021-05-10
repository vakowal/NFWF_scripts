# tabular data for Guam from EPA

public_sec_df <- read.csv("D:/NFWF_PhaseIII/Guam/Data - EPA/bsp_layers93/PUBLIC_SECTOR.csv")
unique(public_sec_df$CAT_SPEC)
cat_spec_summary <- aggregate(ID~CAT_SPEC, data=public_sec_df, FUN=length)
colnames(cat_spec_summary)[2] <- 'number_features'
write.csv(cat_spec_summary, "D:/NFWF_PhaseIII/Guam/Data - EPA/bsp_layers93/PUBLIC_SECTOR_category_summary.csv")

