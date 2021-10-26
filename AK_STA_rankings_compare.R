# Compare rankings of community footprints according to STA combined risk ratings
# and according to our threat and exposure indices
library(ggplot2)
library(GGally)

# rankings according to zonal stats on our indices
nfwf_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v1_revised1_exposure_v2.csv"
nfwf_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v2_exposure_v3.csv"
nfwf_df <- read.csv(nfwf_rank_path)

sta_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Combined_Risk_Ratings_filter.csv"
sta_df <- read.csv(sta_rank_path)
sta_df <- sta_df[, c("ï..FID", "Community", "Combined_r")]
colnames(sta_df) <- c('fid', 'Community', 'STA_combined_risk_rating_rank')

rank_df <- merge(nfwf_df, sta_df)
cor.test(rank_df$STA_combined_risk_rating_rank, rank_df$threat_rank, method='pearson')
cor.test(rank_df$STA_combined_risk_rating_rank, rank_df$exposure_rank, method='pearson')

#parallel coordinates plot
p <- ggparcoord(rank_df, columns=c(8, 5), scale='globalminmax')
p <- p + scale_y_reverse() + ylab("Rank") + xlab("")
print(p)
pngname <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_threat_v2.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()

p <- ggplot(rank_df, aes(x=STA_combined_risk_rating_rank, y=threat_rank))
p <- p + geom_point() + geom_text(x=16, y=125, label="r = 0.61") # 0.76 for threat v2
p <- p + xlab("Rank: STA combined risk rating") + ylab("Rank: Threat Index v1")
pngname <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/scatter_rank_threat_v1_revised1.png" # v1_revised1.png"
png(filename=pngname, width=4, height=4, units='in', res=300)
print(p)
dev.off()

p <- ggplot(rank_df, aes(x=STA_combined_risk_rating_rank, y=threat_mean))
p <- p + geom_point()
print(p)

p <- ggplot(rank_df, aes(x=STA_combined_risk_rating_rank, y=exposure_rank))
p <- p + geom_point() + geom_text(x=16, y=125, label="r = 0.71")
p <- p + xlab("Rank: STA combined risk rating") + ylab("Rank: Exposure Index")
pngname <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/scatter_rank_exposure_v3.png"  # v1_revised1.png"
png(filename=pngname, width=4, height=4, units='in', res=300)
print(p)
dev.off()

p <- ggplot(rank_df, aes(x=STA_combined_risk_rating_rank, y=exposure_mean))
p <- p + geom_point()
print(p)

# rankings of individual threat inputs
nfwf_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_inputs_v2.csv"
nfwf_df <- read.csv(nfwf_rank_path)

sta_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Groups_Combined_Risk_Ratings_filter.csv"
sta_df <- read.csv(sta_rank_path)
sta_df <- sta_df[, c(
    "ï..FID", "Community", "Combined_r", "EROSION_rank", "FLOOD_rank",
    "PERMAFROST_rank", "erosion_group", "flood_group", "permafrost_group")]
colnames(sta_df) <- c(
    'fid', 'Community', 'STA_combined_risk_rating_rank', 'Erosion (STA)',
    'Flood (STA)', 'Permafrost (STA)', 'erosion_group', 'flood_group',
    'permafrost_group')

rank_df <- merge(nfwf_df, sta_df)
colnames(rank_df)[6] <- 'Erosion input (NFWF)'
colnames(rank_df)[7] <- 'Floodprone input (NFWF)'
colnames(rank_df)[8] <- 'Permafrost input (NFWF)'
cor.test(rank_df$sta_erosion_rank, rank_df$erosion_rank, method='pearson')
cor.test(rank_df$sta_flood_rank, rank_df$flooding_rank, method='pearson')
cor.test(rank_df$sta_permafrost_rank, rank_df$permafrost_rank, method='pearson')

rank_df$erosion_group <- factor(rank_df$erosion_group)
rank_df$flood_group <- factor(rank_df$flood_group)
rank_df$permafrost_group <- factor(rank_df$permafrost_group)

# parallel coordinates plots:
p <- ggparcoord(
  rank_df, columns=c(11, 6), groupColumn='erosion_group', scale='globalminmax')
p <- p + scale_y_reverse() + ylab("Rank") + xlab("")
print(p)
pngname <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_erosion_v2.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()

p <- ggparcoord(
  rank_df, columns=c(12, 7), groupColumn='flood_group', scale='globalminmax')
p <- p + scale_y_reverse() + ylab("Rank") + xlab("")
print(p)
pngname <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_flood_v2.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()

p <- ggparcoord(
  rank_df, columns=c(13, 8), groupColumn='permafrost_group', scale='globalminmax')
p <- p + scale_y_reverse() + ylab("Rank") + xlab("")
print(p)
pngname <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_permafrost_v2.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()

# combine individual group and combined risk ratings
sta_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Combined_Risk_Ratings_filter.csv"
sta_df <- read.csv(sta_rank_path)
sta_groups <- read.csv("D:/NFWF_PhaseII/Alaska_archive/Data_Advisory_Committee/State_Threat_Assessment_Groups/STA_Individual_Threat_Rankings.csv")
matched <- merge(sta_df, sta_groups, by.x='Community_', by.y='Community')
write.csv(matched, "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Groups_Combined_Risk_Ratings_filter.csv",
          row.names=FALSE)
