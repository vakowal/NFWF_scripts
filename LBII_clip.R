# crop local biodiversity intactness raster 
library(raster)
library(sf)

GL_boundary <- "D:/NFWF_PhaseIII/Preliminary_Boundaries/GL_PrelimBndy.shp"
AS_boundary <- "D:/NFWF_PhaseIII/Preliminary_Boundaries/AS_PrelimBoundary.shp"
GU_boundary <- "D:/NFWF_PhaseIII/Preliminary_Boundaries/Guam_PrelimBndy.shp"
AK_boundary <- "D:/Packages/AK_Draft_Indices_012721_6145f6/commondata/boundaries/AK_20mDepth_Boundary_v1.shp"

lbii_path <- "D:/Datasets/Local_Biodiversity_Intactness/second_download_6-15-21/lbii.asc"
# fli_path <- "D:/Datasets/Forest_Landscape_Integrity/forest_landscape_integrity_earth.tif"
thepath <- lbii_path
shppath <- GU_boundary  # AS_boundary  # GU boundary, AK boundary
lbii_ras <- raster(thepath)

# set projection of LBII (assume WGS84 lat long)
projection(lbii_ras) <- CRS('+init=EPSG:4326')

# project NFWF boundary to match
pr_boundary <- st_read(shppath)
pr_wgs84 <- st_transform(pr_boundary, crs=CRS('+init=EPSG:4326'))

# crop LBII to NFWF boundary
lbii_pr <- crop(x=lbii_ras, y=pr_wgs84)
# out_path <- "D:/Datasets/Local_Biodiversity_Intactness/lbii_Great_Lakes.tif"  # GU, AS, AK
out_path <- "D:/Datasets/Local_Biodiversity_Intactness/lbii_AK.tif"  # GU, AS, AK
writeRaster(lbii_pr, out_path, format='GTiff')
