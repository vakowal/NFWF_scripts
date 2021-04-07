library(ggplot2)
# try to classify community footprints
fprint_df <- read.csv(
  "D:/Current/Alaska/community_types_exploring/community_footprints_copy.csv")
bldg_count <- read.csv(
  "D:/Current/Alaska/community_types_exploring/community_footprints_building_count.csv")
pop_df <- read.csv(
  "D:/Current/Alaska/community_types_exploring/Alaska_cities_population.csv")
sta_comm_df <- read.csv(
  "D:/packaged_by_Kim_1-14-21/Alaska_1_14_2021/Data_Advisory_Committee/AK_Native_Villages/186 Statewide Threat Assessment Communities.csv")
gao_anv_df <- read.csv(
  "D:/packaged_by_Kim_1-14-21/Alaska_1_14_2021/Data_Advisory_Committee/AK_Native_Villages/GAO 213 ANVs.csv")
stat_areas_df <- read.csv(
  "D:/packaged_by_Kim_1-14-21/Alaska_1_14_2021/Data_Advisory_Committee/AK_Native_Villages/Alaska Native Village Statistical Areas - 2010 Census.csv")

# match the GAO ANV list, which included 213 communities, with census areas
gao_merge <- merge(
  stat_areas_df, gao_anv_df, by.x='BIA.recognized.name.8',
  by.y='Join.field', all.y=TRUE)
write.csv(gao_merge, "D:/packaged_by_Kim_1-14-21/Alaska_1_14_2021/Data_Advisory_Committee/AK_Native_Villages/GAO 213 ANVs_coordinates.csv")

# make sure the 213 GAO communities include all communities in the STA

sta_df <- merge(
  fprint_df, sta_comm_df, by.x='NAMELSAD', by.y='COMMUNITY', all=TRUE)

fprint_df <- merge(fprint_df, bldg_count)
fprint_df$in_STA <- factor(fprint_df$in_STA)
fprint_df$bldg_point_density <- fprint_df$Point_Coun / fprint_df$SHAPE_STAr
p <- ggplot(fprint_df, aes(bldg_point_density, fill=CLASSFP))
p <- p + geom_density(alpha=0.7)
print(p)

p <- ggplot(fprint_df, aes(bldg_point_density, fill=in_STA))
p <- p + geom_density(alpha=0.7)
print(p)

p <- ggplot(fprint_df, aes(x=CLASSFP, y=Point_Coun))
p <- p + geom_boxplot() + ylab("Building footprints point count")
pngname <- "D:/Current/Alaska/community_types_exploring/bldg_footprint_point_count_by_CLASSFP.png"
png(filename=pngname, width=4, height=4, units='in', res=300)
print(p)
dev.off()

community_df <- merge(
  fprint_df, pop_df, by.x='NAMELSAD', by.y='Name', all=TRUE)
community_df[is.na(community_df$in_STA), 'in_STA'] <- 0

no_urb <- community_df[community_df$Population.2010 < 8000, ]
p <- ggplot(community_df, aes(x=CLASSFP, y=Population.2010))
p <- p + geom_boxplot()
print(p)

p <- ggplot(no_urb, aes(Population.2010, fill=in_STA))
p <- p + geom_density(alpha=0.7)
print(p)
pngname <- "D:/Current/Alaska/community_types_exploring/population_by_STA_inclusion_exclude_pop_gt_8000.png"
png(filename=pngname)
print(p)
dev.off()
