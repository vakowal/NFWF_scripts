# Compare rankings of community footprints according to STA combined risk ratings
# and according to our threat and exposure indices
library(ggplot2)
library(GGally)

# rankings according to zonal stats on our indices
nfwf_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v1_revised1_exposure_v2.csv"
nfwf_rank_path <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/zonalmean_threat_v2.csv"
nfwf_rank_path <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v3.csv"
nfwf_df <- read.csv(nfwf_rank_path)

sta_rank_path <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Groups_Combined_Risk_Ratings.csv"
sta_df <- read.csv(sta_rank_path)

rank_df <- merge(nfwf_df, sta_df, by='NAMELSAD')
# calculate ranks from NFWF indices, counting only communities assessed by STA
rank_df$threat_rank <- rank(-rank_df$threat_mean, ties.method='min')
rank_df$erosion_rank <- rank(-rank_df$erosion_mean, ties.method='min')
rank_df$flooding_rank <- rank(-rank_df$flooding_mean, ties.method='min')
rank_df$permafrost_rank <- rank(-rank_df$permafrost_mean, ties.method='min')
cor.test(rank_df$Combined_risk_rating, rank_df$threat_rank, method='pearson')
cor.test(rank_df$Combined_risk_rating, rank_df$exposure_rank, method='pearson')

#parallel coordinates plot
col_idx <- c(
  match('Combined_risk_rating', colnames(rank_df)),
  match('threat_rank', colnames(rank_df)))
p <- ggparcoord(rank_df, columns=col_idx, scale='globalminmax')
p <- p + scale_y_reverse() + ylab("Rank") + xlab("")
print(p)
pngname <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_threat_v3.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()

p <- ggplot(rank_df, aes(x=Combined_risk_rating, y=threat_rank))
p <- p + geom_point() + geom_text(x=16, y=125, label="r = 0.91") # 0.76 for threat v2; 0.61 for v1
p <- p + xlab("Rank: STA combined risk rating") + ylab("Rank: Threat Index v3")
pngname <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/scatter_rank_threat_v4.png"
png(filename=pngname, width=4, height=4, units='in', res=300)
print(p)
dev.off()

p <- ggplot(rank_df, aes(x=Combined_risk_rating, y=threat_mean))
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
# rank_df generated above contains these inputs
# nfwf_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_inputs_v2.csv"
# nfwf_df <- read.csv(nfwf_rank_path)
# 
# sta_rank_path <- "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Groups_Combined_Risk_Ratings_filter.csv"
# sta_df <- read.csv(sta_rank_path)
# sta_df <- sta_df[, c(
#     "ï..FID", "Community", "Combined_r", "EROSION_rank", "FLOOD_rank",
#     "PERMAFROST_rank", "erosion_group", "flood_group", "permafrost_group")]
# colnames(sta_df) <- c(
#     'fid', 'Community', 'STA_combined_risk_rating_rank', 'Erosion (STA)',
#     'Flood (STA)', 'Permafrost (STA)', 'erosion_group', 'flood_group',
#     'permafrost_group')
# 
# rank_df <- merge(nfwf_df, sta_df)
colnames(rank_df)[match('erosion_rank', colnames(rank_df))] <- 'Erosion input (NFWF)'
colnames(rank_df)[match('flooding_rank', colnames(rank_df))] <- 'Floodprone input (NFWF)'
colnames(rank_df)[match('permafrost_rank', colnames(rank_df))] <- 'Permafrost input (NFWF)'
colnames(rank_df)[match('EROSION_rank', colnames(rank_df))] <- 'Erosion (STA)'
colnames(rank_df)[match('FLOOD_rank', colnames(rank_df))] <- 'Flood (STA)'
colnames(rank_df)[match('PERMAFROST_rank', colnames(rank_df))] <- 'Permafrost (STA)'

cor.test(rank_df[, 'Erosion (STA)'], rank_df[, 'Erosion input (NFWF)'], method='pearson')
cor.test(rank_df[, 'Flood (STA)'], rank_df[, 'Floodprone input (NFWF)'], method='pearson')
cor.test(rank_df[, 'Permafrost (STA)'], rank_df[, 'Permafrost input (NFWF)'], method='pearson')
write.csv(rank_df, "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/threat_inputs_ranks.csv",
          row.names=FALSE)

rank_df$erosion_group <- factor(rank_df$erosion_group)
rank_df$flood_group <- factor(rank_df$flood_group)
rank_df$permafrost_group <- factor(rank_df$permafrost_group)

# parallel coordinates plots:
col_idx <- c(
  match('Erosion (STA)', colnames(rank_df)),
  match('Erosion input (NFWF)', colnames(rank_df)))
p <- ggparcoord(
  rank_df, columns=col_idx, groupColumn='erosion_group', scale='globalminmax')
p <- p + scale_y_reverse() + ylab("Rank") + xlab("")
print(p)
pngname <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_erosion_v3.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()

col_idx <- c(
  match('Flood (STA)', colnames(rank_df)),
  match('Floodprone input (NFWF)', colnames(rank_df)))
p <- ggparcoord(
  rank_df, columns=col_idx, groupColumn='flood_group', scale='globalminmax')
p <- p + scale_y_reverse() + ylab("Rank") + xlab("")
print(p)
pngname <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_flood_v3.png"
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
nfwf_rank_path_v1 <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v1_revised1_exposure_v2.csv"
nfwf_rank_path_v3 <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v3.csv"
nfwf_df_v1 <- read.csv(nfwf_rank_path_v1)
# calculate ranks from zonal threat
nfwf_df_v1$threat_rank_v1 <- rank(-nfwf_df_v1$threat_mean, ties.method='min')
nfwf_df_v1 <- nfwf_df_v1[, c('fid', 'threat_rank_v1')]
nfwf_df_v3 <- read.csv(nfwf_rank_path_v3)
nfwf_df_v3$threat_rank_v3 <- rank(-nfwf_df_v3$threat_mean, ties.method='min')
nfwf_df_v3 <- nfwf_df_v3[, c('fid', 'threat_rank_v3')]
nfwf_df <- merge(nfwf_df_v1, nfwf_df_v3)
# link fid to NAMELSAD
fid_match <- read.csv("E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v3.csv")
nfwf_df <- merge(nfwf_df, fid_match) 
sta_rank_path <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Groups_Combined_Risk_Ratings.csv"
sta_df <- read.csv(sta_rank_path)
sta_df <- sta_df[, c('NAMELSAD', 'Combined_risk_rating')]
rank_df <- merge(nfwf_df, sta_df, by='NAMELSAD')
write.csv(rank_df, "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/STA_threat_v1_v3_rankings.csv",
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

# assign threat input ranks in {1:5}: permafrost
permaf_gr1_split <- quantile(
  fullmatch_df[fullmatch_df$permafrost_group == 1, 'PERMAFROST_rank'],
  probs=0.5, na.rm=TRUE)[1]
permaf_gr2_split <- quantile(
  fullmatch_df[fullmatch_df$permafrost_group == 2, 'PERMAFROST_rank'],
  probs=0.5, na.rm=TRUE)[1]

fullmatch_df[
  !(is.na(fullmatch_df$permafrost_group)) &
    (fullmatch_df$PERMAFROST_rank < permaf_gr1_split), 'pf_in_rank'] <- 5
fullmatch_df[
  !(is.na(fullmatch_df$permafrost_group)) &
    (fullmatch_df$permafrost_group == 1) &
    (fullmatch_df$PERMAFROST_rank >= permaf_gr1_split), 'pf_in_rank'] <- 4
fullmatch_df[
  !(is.na(fullmatch_df$permafrost_group)) &
    (fullmatch_df$permafrost_group == 2) & 
    (fullmatch_df$PERMAFROST_rank < permaf_gr2_split), 'pf_in_rank'] <- 3
fullmatch_df[
  !(is.na(fullmatch_df$permafrost_group)) &
    (fullmatch_df$permafrost_group == 2) &
    (fullmatch_df$PERMAFROST_rank >= permaf_gr2_split), 'pf_in_rank'] <- 2
fullmatch_df[
  !(is.na(fullmatch_df$permafrost_group)) &
    (fullmatch_df$permafrost_group == 3), 'pf_in_rank'] <- 1

in_rank_df <- fullmatch_df[, c("NAMELSAD", "er_in_rank", "fl_in_rank", "pf_in_rank")]
write.csv(
  in_rank_df,
  "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_flood_erosion_permafrost_5group.csv",
  row.names=FALSE)

# make parallel coordinates plot showing the derivation of threat inputs
in_rank_df <- read.csv("E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_flood_erosion_permafrost_5group.csv")
sta_rank_path <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/Community_Footprints_STA_Groups_Combined_Risk_Ratings.csv"
sta_df <- read.csv(sta_rank_path)
merged_df <- merge(in_rank_df, sta_df, by='NAMELSAD')
colnames(merged_df)[1:9] <- c(
  'NAMELSAD', 'erosion input', 'floodprone input', 'permafrost input',
  'Community',  'combined risk rating', 'erosion rank (STA)',
  'flood rank (STA)', 'permafrost rank (STA)')

merged_df$erosion_group <- factor(merged_df$erosion_group)
merged_df$flood_group <- factor(merged_df$flood_group)
merged_df$permafrost_group <- factor(merged_df$permafrost_group)

col_idx <- c(
  match('erosion rank (STA)', colnames(merged_df)),
  match('erosion input', colnames(merged_df)))
p <- ggparcoord(merged_df, columns=col_idx, groupColumn='erosion_group', scale='std')
p <- p + scale_y_reverse() + ylab("Normalized rank") + xlab("")
print(p)
pngname <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_erosion_input_v3.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()

col_idx <- c(
  match('flood rank (STA)', colnames(merged_df)),
  match('floodprone input', colnames(merged_df)))
p <- ggparcoord(merged_df, columns=col_idx, groupColumn='flood_group', scale='std')
p <- p + scale_y_reverse() + ylab("Normalized rank") + xlab("")
print(p)
pngname <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_flood_input_v3.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()

col_idx <- c(
  match('permafrost rank (STA)', colnames(merged_df)),
  match('permafrost input', colnames(merged_df)))
p <- ggparcoord(merged_df, columns=col_idx, groupColumn='permafrost_group', scale='std')
p <- p + scale_y_reverse() + ylab("Normalized rank") + xlab("")
print(p)
pngname <- "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/parcoord_permafrost_input_v3.png"
png(filename=pngname, width=5, height=4, units='in', res=300)
print(p)
dev.off()
