"""Shapefile operations for Great Lakes NFWF indices."""
import os

import geopandas

# NCRF boundary for Great Lakes
_BOUNDARY_PATH = "D:/NFWF_PhaseIII/Great_Lakes/GreatLakes_NRCFBoundary/NCRF_Coastal_Areas_GreatLakes.shp"


def process_NWI():
    """Process NWI shapefiles across states."""
    NWI_path_list = [
        "D:/Datasets/FWS/NWI/Great_Lakes/WI_shapefile_wetlands/WI_Wetlands_South.shp",
        "D:/Datasets/FWS/NWI/Great_Lakes/WI_shapefile_wetlands/WI_Wetlands_North.shp",
        "D:/Datasets/FWS/NWI/Great_Lakes/PA_shapefile_wetlands/PA_shapefile_wetlands/PA_Wetlands.shp",
        "D:/Datasets/FWS/NWI/Great_Lakes/OH_shapefile_wetlands/OH_shapefile_wetlands/OH_Wetlands.shp",
        "D:/Datasets/FWS/NWI/Great_Lakes/NY_shapefile_wetlands/NY_Wetlands.shp",
        "D:/Datasets/FWS/NWI/Great_Lakes/MI_shapefile_wetlands/MI_shapefile_wetlands/MI_Wetlands.shp",
        "D:/Datasets/FWS/NWI/Great_Lakes/IN_shapefile_wetlands/IN_Wetlands.shp",
        "D:/Datasets/FWS/NWI/Great_Lakes/IL_shapefile_wetlands/IL_Wetlands.shp"]

    # project the boundary to match the NWI data
    example_NWI_gdf = geopandas.read_file(NWI_path_list[0])
    target_crs = example_NWI_gdf.crs
    boundary_gdf = geopandas.read_file(_BOUNDARY_PATH)
    boundary_proj = boundary_gdf.to_crs(target_crs)

    # clip each NWI dataset to boundary
    # ran into problem here:
    # TopologyException: Input geom 0 is invalid: Ring Self-intersection at or near point
    merged_gdf = []
    for NWI_path in NWI_path_list:
        NWI_gdf = geopandas.read_file(NWI_path)
        # NWI_gdf.sindex.query(boundary_proj)  ?? OR
        NWI_clip = geopandas.clip(NWI_gdf, boundary_proj)
        merged_gdf.append(NWI_clip)

    # write merged NWI across states
    merged_path = "D:/NFWF_PhaseIII/Great_Lakes/Data - wildlife/NWI/NWI_merged.shp"
    merged_gdf.to_file(merged_path)


if __name__ == "__main__":
    process_NWI()