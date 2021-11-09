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
cor.test(rank_df[, 'Erosion (STA)'], rank_df[, 'Erosion input (NFWF)'], method='pearson')
cor.test(rank_df[, 'Flood (STA)'], rank_df[, 'Floodprone input (NFWF)'], method='pearson')
cor.test(rank_df[, 'Permafrost (STA)'], rank_df[, 'Permafrost input (NFWF)'], method='pearson')
write.csv(rank_df, "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/threat_inputs_ranks.csv",
          row.names=FALSE)

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

# make table of rankings by community
nfwf_rank_path_v1 <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v1_revised1_exposure_v2.csv"
nfwf_rank_path_v2 <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v2_exposure_v3.csv"
nfwf_df_v1 <- read.csv(nfwf_rank_path_v1)
nfwf_df_v1 <- nfwf_df_v1[, c('fid', 'threat_rank', 'exposure_rank')]
colnames(nfwf_df_v1) <- c('fid', 'threat_rank_v1', 'exposure_rank_v1')
nfwf_df_v2 <- read.csv(nfwf_rank_path_v2)
nfwf_df_v2 <- nfwf_df_v2[, c('fid', 'threat_rank', 'exposure_rank')]
colnames(nfwf_df_v2) <- c('fid', 'threat_rank_v2', 'exposure_rank_v2')
nfwf_df <- merge(nfwf_df_v1, nfwf_df_v2)
sta_df <- read.csv("D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Combined_Risk_Ratings_filter.csv")
sta_df <- sta_df[, c("ï..FID", "Community_", "Combined_r")]
colnames(sta_df) <- c('fid', 'Community', 'STA_combined_risk_rating_rank')
rank_df <- merge(nfwf_df, sta_df)
write.csv(rank_df, "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/STA_threat_exposure_v1_v2_rankings.csv",
          row.names=FALSE)

# generate scores from the STA to be substituted for erosion and flooding inputs to the threat index
# assign each community a value in 1-5 according to its rank from the STA.
# This table contains combined rankings, individual rankings, and STA-assigned groups for all communities
# it also contains a field "NAMELSAD" which can be used to join it to spatial data for community footprints
fullmatch <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Groups_Combined_Risk_Ratings.csv"
fullmatch_df <- read.csv(fullmatch)

# assign threat input ranks in {1:5}: erosion
erosion_gr1_split <- quantile(
  fullmatch_df[fullmatch_df$erosion_group == 1, 'EROSION_rank'],
  probs=0.5, na.rm=TRUE)[1]
erosion_gr2_split <- quantile(
  fullmatch_df[fullmatch_df$erosion_group == 2, 'EROSION_rank'],
  probs=0.5, na.rm=TRUE)[1]

fullmatch_df[
  !(is.na(fullmatch_df$erosion_group)) &
  (fullmatch_df$EROSION_rank < erosion_gr1_split), 'er_in_rank'] <- 5
fullmatch_df[
  !(is.na(fullmatch_df$erosion_group)) &
  (fullmatch_df$erosion_group == 1) &
  (fullmatch_df$EROSION_rank >= erosion_gr1_split), 'er_in_rank'] <- 4
fullmatch_df[
  !(is.na(fullmatch_df$erosion_group)) &
  (fullmatch_df$erosion_group == 2) & 
  (fullmatch_df$EROSION_rank < erosion_gr2_split), 'er_in_rank'] <- 3
fullmatch_df[
  !(is.na(fullmatch_df$erosion_group)) &
  (fullmatch_df$erosion_group == 2) &
  (fullmatch_df$EROSION_rank >= erosion_gr2_split), 'er_in_rank'] <- 2
fullmatch_df[
  !(is.na(fullmatch_df$erosion_group)) &
  (fullmatch_df$erosion_group == 3), 'er_in_rank'] <- 1

# assign threat input ranks in {1:5}: flooding
flooding_gr1_split <- quantile(
  fullmatch_df[fullmatch_df$flood_group == 1, 'FLOOD_rank'],
  probs=0.5, na.rm=TRUE)[1]
flooding_gr2_split <- quantile(
  fullmatch_df[fullmatch_df$flood_group == 2, 'FLOOD_rank'],
  probs=0.5, na.rm=TRUE)[1]

fullmatch_df[
  !(is.na(fullmatch_df$flood_group)) &
  (fullmatch_df$FLOOD_rank < flooding_gr1_split), 'fl_in_rank'] <- 5
fullmatch_df[
  !(is.na(fullmatch_df$flood_group)) &
  (fullmatch_df$flood_group == 1) &
    (fullmatch_df$FLOOD_rank >= flooding_gr1_split), 'fl_in_rank'] <- 4
fullmatch_df[
  !(is.na(fullmatch_df$flood_group)) &
  (fullmatch_df$flood_group == 2) & 
    (fullmatch_df$FLOOD_rank < flooding_gr2_split), 'fl_in_rank'] <- 3
fullmatch_df[
  !(is.na(fullmatch_df$flood_group)) &
  (fullmatch_df$flood_group == 2) &
    (fullmatch_df$FLOOD_rank >= flooding_gr2_split), 'fl_in_rank'] <- 2
fullmatch_df[
  !(is.na(fullmatch_df$flood_group)) &
  (fullmatch_df$flood_group == 3), 'fl_in_rank'] <- 1
in_rank_df <- fullmatch_df[, c("NAMELSAD", "er_in_rank", "fl_in_rank")]
write.csv(
  in_rank_df,
  "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_flood_erosion_5group.csv",
  row.names=FALSE)
