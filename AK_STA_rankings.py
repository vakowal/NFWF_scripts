"""Compare community rankings from NFWF indices to STA."""
import os

import pandas

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing

# threat index path
_THREAT_PATH = "D:/Packages/AK_Threat_Inputs_012721_9522b9/AK_Threat_Index_10class_revised1.tif"

# exposure index path
_EXPOSURE_PATH = "D:/Packages/AK_Asset_Inputs_012721_2fd5ef/commondata/raster_data/AK_Exposure_Index_10class_v2.tif"

# community footprints path
_FOOTPRINTS_PATH = "D:/NFWF_PhaseII/Alaska/community_types_exploring/Community_Footprints_STA_Combined_Risk_Ratings_filter.shp"


def calc_rankings(save_as):
    """Calculate zonal mean index inside community footprints.

    Args:
        save_as (string): the path to save zonal mean statistics.

    """
    zonal_dict = {
        'fid': [],
        'threat_mean': [],
        'exposure_mean': [],
    }
    
    # zonal mean of threat
    threat_stats = pygeoprocessing.zonal_statistics(
        (_THREAT_PATH, 1), _FOOTPRINTS_PATH)
    # zonal mean of exposure
    exposure_stats = pygeoprocessing.zonal_statistics(
        (_EXPOSURE_PATH, 1), _FOOTPRINTS_PATH)
    
    for fid in threat_stats:
        zonal_dict['fid'].append(fid)
        zonal_dict['threat_mean'].append(
            (threat_stats[fid]['sum'] / threat_stats[fid]['count']))
        zonal_dict['exposure_mean'].append(
            (exposure_stats[fid]['sum'] / exposure_stats[fid]['count']))
    zonal_df = pandas.DataFrame.from_dict(zonal_dict, orient='columns')
    zonal_df['threat_rank'] = zonal_df['threat_mean'].rank(ascending=False)
    zonal_df['exposure_rank'] = zonal_df['exposure_mean'].rank(ascending=False)
    zonal_df.to_csv(save_as)


if __name__ == "__main__":
    zonal_table_path = "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat v1_revised1_exposure_v2.csv"
    calc_rankings(zonal_table_path)