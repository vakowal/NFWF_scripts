# crop local biodiversity intactness raster 
library(raster)
library(sf)

thepath <- "D:/Datasets/Local_Biodiversity_Intactness/lbii.asc"
shppath <- "D:/NFWF_PhaseIII/Preliminary_Boundaries/GL_PrelimBndy.shp"
lbii_ras <- raster(thepath)

# set projection of LBII (assume WGS84 lat long)
projection(lbii_ras) <- CRS('+init=EPSG:4326')

# project PR boundary to match
pr_boundary <- st_read(shppath)
pr_wgs84 <- st_transform(pr_boundary, crs=CRS('+init=EPSG:4326'))

# crop LBII to PR boundary
lbii_pr <- crop(x=lbii_ras, y=pr_wgs84)
out_path <- "D:/Datasets/Local_Biodiversity_Intactness/lbii_Great_Lakes.tif"
writeRaster(lbii_pr, out_path, format='GTiff')
