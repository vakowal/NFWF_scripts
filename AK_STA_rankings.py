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
# _THREAT_PATH = "E:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_Threat_Index_v2.tif"
# v3
# _THREAT_PATH = "E:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v3_110921/AK_Threat_Index_v3.tif"
# v4
_THREAT_PATH = "E:/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data11/AK_Threat_Index_v4.tif"

# exposure index path
# v2
# _EXPOSURE_PATH = "D:/Packages/AK_Asset_Inputs_012721_2fd5ef/commondata/raster_data/AK_Exposure_Index_10class_v2.tif"
# v3 (calculated from threat index v2)
# _EXPOSURE_PATH = "D:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_Exposure_Index_v3.tif"
# v5 (calculated from threat index v4 and asset index v3)
_EXPOSURE_PATH = "E:/Packages/AK_Asset_Inputs_012721_2fd5ef/commondata/raster_data1/AK_Exposure_Index_10class_v5.tif"

# community footprints path
_FOOTPRINTS_PATH = "E:/NFWF_PhaseII/Alaska/community_types_exploring/Community_Footprints_STA_Groups_Combined_Risk_Ratings_aoi.shp"


def calc_rankings(save_as):
    """Calculate zonal mean index inside community footprints.

    Args:
        save_as (string): the path to save zonal mean statistics.

    """
    zonal_dict = {
        'fid': [],
        # 'threat_mean': [],
        'exposure_mean': [],
    }
    
    # zonal mean of threat
    # print("calculating zonal threat")
    # threat_stats = pygeoprocessing.zonal_statistics(
    #     (_THREAT_PATH, 1), _FOOTPRINTS_PATH)
    # zonal mean of exposure
    print("calculating zonal exposure")
    exposure_stats = pygeoprocessing.zonal_statistics(
        (_EXPOSURE_PATH, 1), _FOOTPRINTS_PATH)
    
    for fid in exposure_stats:
        try:
            exposure_mean = (
                exposure_stats[fid]['sum'] / exposure_stats[fid]['count'])
        except ZeroDivisionError:
            # do not include communities that lie outside threat index
            continue
        zonal_dict['fid'].append(fid)
        # zonal_dict['threat_mean'].append(
        #     (threat_stats[fid]['sum'] / threat_stats[fid]['count']))
        zonal_dict['exposure_mean'].append(exposure_mean)
    zonal_df = pandas.DataFrame.from_dict(zonal_dict, orient='columns')

    # add NAMELSAD field so that these can be compared to STA table
    aggregate_vector = gdal.OpenEx(_FOOTPRINTS_PATH, gdal.OF_VECTOR)
    aggregate_layer = aggregate_vector.GetLayer()
    fid_list = [feature.GetFID() for feature in aggregate_layer]
    name_list = [feature.GetField('NAMELSAD') for feature in aggregate_layer]
    match_dict = {'fid': fid_list, 'NAMELSAD': name_list}
    match_df = pandas.DataFrame.from_dict(match_dict, orient='columns')
    zonal_plus_fid = zonal_df.merge(
        match_df, how='outer', on='fid', suffixes=(None, None))
    zonal_plus_fid.to_csv(save_as)


def input_rankings(save_as):
    """Calculate zonal mean of threat index inputs.

    Args:
        save_as (string): the path to save zonal mean statistics.

    """
    # erosion_path = "D:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_erosion_v2_rc.tif"
    # flooding_path = "D:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_floodprone_v2.tif"
    erosion_path = "E:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v3_110921/AK_erosion_v3.tif"
    flooding_path = "E:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v3_110921/AK_floodprone_v3.tif"
    permafrost_path = "E:/NFWF_PhaseII/Alaska/Revise_threat_index/threat_v2_101421/AK_permafrost_v2.tif"
    zonal_dict = {
        'fid': [],
        'threat_mean': [],
        # 'erosion_mean': [],
        # 'flooding_mean': [],
        # 'permafrost_mean': [],
    }
    
    # zonal mean of threat
    print("calculating zonal threat")
    threat_stats = pygeoprocessing.zonal_statistics(
        (_THREAT_PATH, 1), _FOOTPRINTS_PATH)
    # zonal mean of erosion
    # print("calculating zonal erosion")
    # erosion_stats = pygeoprocessing.zonal_statistics(
    #     (erosion_path, 1), _FOOTPRINTS_PATH)
    # # zonal mean of flooding
    # print("calculating zonal flooding")
    # flooding_stats = pygeoprocessing.zonal_statistics(
    #     (flooding_path, 1), _FOOTPRINTS_PATH)
    # # zonal mean permafrost
    # print("calculating zonal permafrost")
    # permafrost_stats = pygeoprocessing.zonal_statistics(
    #     (permafrost_path, 1), _FOOTPRINTS_PATH)
    
    for fid in threat_stats:
        try:
            threat_mean = threat_stats[fid]['sum'] / threat_stats[fid]['count']            
        except ZeroDivisionError:
            # do not include communities that lie outside threat index
            continue
        zonal_dict['fid'].append(fid)
        zonal_dict['threat_mean'].append(threat_mean)
        # zonal_dict['erosion_mean'].append(
        #     (erosion_stats[fid]['sum'] / erosion_stats[fid]['count']))
        # zonal_dict['flooding_mean'].append(
        #     (flooding_stats[fid]['sum'] / flooding_stats[fid]['count']))
        # zonal_dict['permafrost_mean'].append(
        #     (permafrost_stats[fid]['sum'] /
        #         permafrost_stats[fid]['count']))
    zonal_df = pandas.DataFrame.from_dict(zonal_dict, orient='columns')
    # zonal_df['threat_rank'] = zonal_df[
    #     'threat_mean'].rank(method='dense', ascending=False)
    # zonal_df['erosion_rank'] = zonal_df[
    #     'erosion_mean'].rank(method='dense', ascending=False)
    # zonal_df['flooding_rank'] = zonal_df[
    #     'flooding_mean'].rank(method='dense', ascending=False)
    # zonal_df['permafrost_rank'] = zonal_df[
    #     'permafrost_mean'].rank(method='dense', ascending=False)

    # add NAMELSAD field so that these can be compared to STA table
    aggregate_vector = gdal.OpenEx(_FOOTPRINTS_PATH, gdal.OF_VECTOR)
    aggregate_layer = aggregate_vector.GetLayer()
    fid_list = [feature.GetFID() for feature in aggregate_layer]
    name_list = [feature.GetField('NAMELSAD') for feature in aggregate_layer]
    match_dict = {'fid': fid_list, 'NAMELSAD': name_list}
    match_df = pandas.DataFrame.from_dict(match_dict, orient='columns')
    zonal_plus_fid = zonal_df.merge(
        match_df, how='outer', on='fid', suffixes=(None, None))
    zonal_plus_fid.to_csv(save_as)


def fid_to_name_table(vector_layer, save_as):
    """Generate table linking FID to NAMELSAD for features in `vector_layer`."""
    aggregate_vector = gdal.OpenEx(vector_layer, gdal.OF_VECTOR)
    aggregate_layer = aggregate_vector.GetLayer()
    fid_list = [feature.GetFID() for feature in aggregate_layer]
    name_list = [feature.GetField('NAMELSAD') for feature in aggregate_layer]
    match_dict = {'fid': fid_list, 'NAMELSAD': name_list}
    match_df = pandas.DataFrame.from_dict(match_dict, orient='columns')
    match_df.to_csv(save_as)


if __name__ == "__main__":
    # zonal_table_path = "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/zonalmean_threat_v2.csv"
    # input_rankings(zonal_table_path)
    # fid_to_name_table(zonal_table_path)
    input_zonal_path = "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/zonalmean_threat_v4.csv"
    # input_rankings(input_zonal_path)
    # fid_to_name_table()
    filtered_footprints = "E:/NFWF_PhaseII/Alaska/community_types_exploring/Community_Footprints_STA_Combined_Risk_Ratings_filter.shp"
    fid_match = "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/fid_to_namelsad_v1.csv"
    fid_to_name_table(filtered_footprints, fid_match)
    zonal_stat_path = "E:/NFWF_PhaseII/Alaska/community_types_exploring/threat_exposure_rankings/zonalmean_exposure_v5.csv"
    # calc_rankings(zonal_stat_path)
