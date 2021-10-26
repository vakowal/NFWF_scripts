"""Compare community rankings from NFWF indices to STA."""
import os

import pandas

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing

# threat index path
# v1_revised_1
# _THREAT_PATH = "D:/Packages/AK_Threat_Inputs_012721_9522b9/AK_Threat_Index_10class_revised1.tif"
# v2
_THREAT_PATH = "D:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_Threat_Index_v2.tif"

# exposure index path
# v2
# _EXPOSURE_PATH = "D:/Packages/AK_Asset_Inputs_012721_2fd5ef/commondata/raster_data/AK_Exposure_Index_10class_v2.tif"
# v3 (calculated from threat index v2)
_EXPOSURE_PATH = "D:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_Exposure_Index_v3.tif"

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
    print("calculating zonal threat")
    threat_stats = pygeoprocessing.zonal_statistics(
        (_THREAT_PATH, 1), _FOOTPRINTS_PATH)
    # zonal mean of exposure
    print("calculating zonal exposure")
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


def input_rankings(save_as):
    """Calculate zonal mean of threat index inputs.

    Args:
        save_as (string): the path to save zonal mean statistics.

    """
    erosion_path = "D:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_erosion_v2_rc.tif"
    flooding_path = "D:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_floodprone_v2.tif"
    permafrost_path = "D:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_permafrost_v2.tif"
    zonal_dict = {
        'fid': [],
        'erosion_mean': [],
        'flooding_mean': [],
        'permafrost_mean': [],
    }
    
    # zonal mean of erosion
    print("calculating zonal erosion")
    erosion_stats = pygeoprocessing.zonal_statistics(
        (erosion_path, 1), _FOOTPRINTS_PATH)
    # zonal mean of flooding
    print("calculating zonal flooding")
    flooding_stats = pygeoprocessing.zonal_statistics(
        (flooding_path, 1), _FOOTPRINTS_PATH)
    # zonal mean permafrost
    print("calculating zonal permafrost")
    permafrost_stats = pygeoprocessing.zonal_statistics(
        (permafrost_path, 1), _FOOTPRINTS_PATH)
    
    for fid in erosion_stats:
        zonal_dict['fid'].append(fid)
        zonal_dict['erosion_mean'].append(
            (erosion_stats[fid]['sum'] / erosion_stats[fid]['count']))
        zonal_dict['flooding_mean'].append(
            (flooding_stats[fid]['sum'] / flooding_stats[fid]['count']))
        zonal_dict['permafrost_mean'].append(
            (permafrost_stats[fid]['sum'] / permafrost_stats[fid]['count']))
    zonal_df = pandas.DataFrame.from_dict(zonal_dict, orient='columns')
    zonal_df['erosion_rank'] = zonal_df['erosion_mean'].rank(ascending=False)
    zonal_df['flooding_rank'] = zonal_df['flooding_mean'].rank(ascending=False)
    zonal_df['permafrost_rank'] = zonal_df[
        'permafrost_mean'].rank(ascending=False)
    zonal_df.to_csv(save_as)


if __name__ == "__main__":
    # zonal_table_path = "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_v2_exposure_v3.csv"
    # calc_rankings(zonal_table_path)
    input_zonal_path = "D:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/rankings_threat_inputs_v2.csv"
    input_rankings(input_zonal_path)