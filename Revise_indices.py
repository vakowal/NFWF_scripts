"""Revise NFWF indices."""
import os
import tempfile
import shutil

import numpy

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing

# items to be added together should have this nodata value
_TARGET_NODATA = 255


def raster_list_sum(
        raster_list, input_nodata, target_path, target_nodata, target_datatype,
        nodata_remove=False):
    """Calculate the sum per pixel across rasters in a list.

    Sum the rasters in `raster_list` element-wise, allowing nodata values
    in the rasters to propagate to the result or treating nodata as zero. If
    nodata is treated as zero, areas where all inputs are nodata will be nodata
    in the output.

    Parameters:
        raster_list (list): list of paths to rasters to sum
        input_nodata (float or int): nodata value in the input rasters
        target_path (string): path to location to store the result
        target_nodata (float or int): nodata value for the result raster
        target_datatype (gdal data type enum): GDAL data type for result raster
        nodata_remove (bool): if true, treat nodata values in input
            rasters as zero. If false, the sum in a pixel where any input
            raster is nodata is nodata.

    Side effects:
        modifies or creates the raster indicated by `target_path`

    Returns:
        None

    """
    def raster_sum_op(*raster_list):
        """Add the rasters in raster_list without removing nodata values."""
        invalid_mask = numpy.any(
            numpy.isclose(numpy.array(raster_list), input_nodata), axis=0)
        for r in raster_list:
            numpy.place(r, numpy.isclose(r, input_nodata), [0])
        sum_of_rasters = numpy.sum(raster_list, axis=0)
        sum_of_rasters[invalid_mask] = target_nodata
        return sum_of_rasters

    def raster_sum_op_nodata_remove(*raster_list):
        """Add the rasters in raster_list, treating nodata as zero."""
        invalid_mask = numpy.all(
            numpy.isclose(numpy.array(raster_list), input_nodata), axis=0)
        for r in raster_list:
            numpy.place(r, numpy.isclose(r, input_nodata), [0])
        sum_of_rasters = numpy.sum(raster_list, axis=0)
        sum_of_rasters[invalid_mask] = target_nodata
        return sum_of_rasters

    if nodata_remove:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in raster_list], raster_sum_op_nodata_remove,
            target_path, target_datatype, target_nodata)

    else:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in raster_list], raster_sum_op,
            target_path, target_datatype, target_nodata)


def revise_index_with_additions(
        new_input_dict, intermediate_dir, boundary_path, existing_index_path,
        revised_index_path):
    """Revise an existing index by adding vector-based features.

    Parameters:
        new_input_dict (dict): a nested python dictionary containing info about 
            new inputs to be added to the index. The keys of the index should
            be the labels for each new input; values in the dictionary should
            be dictionaries containing the keys 'path' (path to the shapefile
            containing new features) and 'rank' (value that each feature should
            be assigned). The shapefile for each input must have the same
            projection as the existing index
        intermediate_dir (string): path to directory where intermediate results
            should be written
        boundary_path (string): path to shapefile giving the boundary of the
            study area, to which the revised index should be clipped
        existing_index_path (string): path to existing index to which new
            additions should be made
        revised_index_path (string): path to location where revised index
            should be saved

    Side effects:
        creates a raster in `intermediate_dir` for each new input
        creates a raster at the location `revised_index_path`

    Returns:
        None

    """
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)
    raster_info = pygeoprocessing.get_raster_info(existing_index_path)
    input_datatype = raster_info['datatype']
    input_nodata = raster_info['nodata'][0]
    pixel_size = raster_info['pixel_size']
    destination_proj = osr.SpatialReference()
    destination_proj.ImportFromWkt(raster_info['projection_wkt'])    

    # ensure that each new input is projected to match existing index
    problem_list = []
    for new_input in new_input_dict:
        new_proj = osr.SpatialReference()
        new_proj.ImportFromWkt(
            pygeoprocessing.get_vector_info(
                new_input_dict[new_input]['path'])['projection_wkt'])
        if (new_proj.IsSame(destination_proj) == 0):
            problem_list.append(new_input)
    if problem_list:
        raise ValueError(
            "Project these to match the existing index: {}".format(
                problem_list))

    # rasterize each new input, using the existing index as a template
    intermediate_path_dict = {}
    for new_input in new_input_dict:
        intermediate_path_dict[new_input] = os.path.join(
            intermediate_dir, '{}.tif'.format(new_input))
        pygeoprocessing.new_raster_from_base(
            existing_index_path, intermediate_path_dict[new_input],
            input_datatype, [input_nodata], fill_value_list=[input_nodata])
        pygeoprocessing.rasterize(
            new_input_dict[new_input]['path'],
            intermediate_path_dict[new_input],
            burn_values=[new_input_dict[new_input]['rank']])

    # add new inputs to existing index; nodata reclassified to zero
    with tempfile.NamedTemporaryFile(
            prefix='revised_index_not_clipped', suffix='.tif',
            delete=False) as unclipped_raster_file:
        unclipped_raster_path = unclipped_raster_file.name
    sum_path_list = [existing_index_path] + [
        intermediate_path_dict[i] for i in intermediate_path_dict]
    raster_list_sum(
        sum_path_list, input_nodata, unclipped_raster_path, input_nodata,
        input_datatype, nodata_remove=True)
    # clip their sum to the region
    pygeoprocessing.align_and_resize_raster_stack(
        [unclipped_raster_path], [revised_index_path], ['near'],
        pixel_size, 'intersection', base_vector_path_list=[boundary_path],
        raster_align_index=0)

    # clean up
    os.remove(unclipped_raster_path)

    
def AS_threat_index_v2():
    """Add tsunami inundation to AS threat index."""
    new_input_dict = {
        'tsunami': {
            'path': "D:/NFWF_PhaseIII/American_Samoa/Data - Exposure Index/Threat_Index/Tsunami/AmericanSamoa_PMTInundationZone_slowWalk_110_timeMap_proj.shp",
            'rank': 3,
        },
    }
    intermediate_dir = "D:/NFWF_PhaseIII/American_Samoa/Data - Exposure Index/Threat_Index/Tsunami"
    boundary_path = "D:/Packages/AS_Threat_Inputs_v1_01252021_de5873/commondata/boundaries/AS_30mDepth_Bndy.shp"
    existing_index_path = "D:/Packages/AS_Threat_Inputs_v1_01252021_de5873/commondata/raster_data10/AS_Threat_Index_v1.tif"
    revised_index_path = "D:/Packages/AS_Threat_Inputs_v1_01252021_de5873/commondata/raster_data10/AS_Threat_Index_v2.tif"
    revise_index_with_additions(
        new_input_dict, intermediate_dir, boundary_path, existing_index_path,
        revised_index_path)


def AS_terrestrial_index_revisions():
    """Add important veg data from Adam Miles to AS terrestrial index."""
    new_input_dict = {
        'primary_forest': {
            'path': "D:/NFWF_PhaseIII/American_Samoa/Data - Wildlife Index_/Wildlife_Index/Data_Adam_Miles/intermediate/primary_forest_merge.shp",
            'rank': 2,
        },
        'naumati_forest': {
            'path': "D:/NFWF_PhaseIII/American_Samoa/Data - Wildlife Index_/Wildlife_Index/Data_Adam_Miles/intermediate/Naumati_forest.shp",
            'rank': 2,
        }
    }
    intermediate_dir = "D:/NFWF_PhaseIII/American_Samoa/Data - Wildlife Index_/Wildlife_Index/Data_Adam_Miles/processed"
    boundary_path = "D:/Packages/AS_Threat_Inputs_v1_01252021_de5873/commondata/boundaries/AS_30mDepth_Bndy.shp"
    existing_index_path = "D:/Packages/AS_Wildlife_Inputs_v1_01252021_52abd4/commondata/raster_data6/AS_Terrestrial_Index_ALL_v1.tif"
    revised_index_path = "D:/Packages/AS_Wildlife_Inputs_v1_01252021_52abd4/commondata/raster_data6/AS_Terrestrial_Index_ALL_v2.tif"
    revise_index_with_additions(
        new_input_dict, intermediate_dir, boundary_path, existing_index_path,
        revised_index_path)


def GU_marine_index_revisions():
    """Add mangrove data to Guam marine index v1."""
    new_input_dict = {
        'mangrove_UnivGuam': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - Wildlife Index/Marine_Index/Mangrove/Univ_Guam_Marine_Lab/mangrove_poly_mg_hfa_proj.shp",
            'rank': 5,
        },
    }
    intermediate_dir = TODO
    boundary_path = "D:/Packages/GU_Wildlife_Inputs_v1_01252021_c4e628/commondata/boundaries/GU_30m_Boundary.shp"
    existing_index_path = "D:/Packages/GU_Wildlife_Inputs_v1_01252021_c4e628/commondata/raster_data1/GU_Marine_Index_All_v1.tif"
    revised_index_path = TODO
    revise_index_with_additions(
        new_input_dict, intermediate_dir, boundary_path, existing_index_path,
        revised_index_path)


def GU_terrestrial_index_revisions():
    """Add data to Guam terrestrial index v1revised1."""
    new_input_dict = {
        'kingfisher_hab': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - Wildlife Index/Kingfisher/Guam_Kingfisher_Potential_Recovery_Habitat_April2015_proj.shp",
            'rank': 2,
        },
        'conservation_areas': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - digital atlas of Southern Guam/protected_areas/infrastructure-protected-areas/conservation_areas_proj.shp",
            'rank': 2,
        },
        'ecological_reserve_areas': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - digital atlas of Southern Guam/protected_areas/infrastructure-protected-areas/ecological_reserve_area_proj.shp",
            'rank': 2,
        }
    }
    intermediate_dir = "D:/NFWF_PhaseIII/Guam/Data - Wildlife Index/terrestrial_index_revisions"
    boundary_path = "D:/Packages/GU_Wildlife_Inputs_v1_01252021_c4e628/commondata/boundaries/GU_30m_Boundary.shp"
    existing_index_path = "D:/Packages/GU_Wildlife_Index_v1revised1_e193be/GU_Terrestrial_Index_ALL_v1revised1.tif"
    revised_index_path = "D:/Packages/GU_Wildlife_Index_v1revised1_e193be/GU_Terrestrial_Index_ALL_v2.tif"
    revise_index_with_additions(
        new_input_dict, intermediate_dir, boundary_path, existing_index_path,
        revised_index_path)


def GU_threat_index_revisions():
    """Add data to Guam threat index v1."""
    new_input_dict = {
        'tsunami': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - Exposure Index/NOAA_tsunami_aware/tsunami_proj.shp",
            'rank': 5,
        },
        'landslides_low': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - Exposure Index/NOAA_tsunami_aware/landslides_low.shp",
            'rank': 1,
        },
        'landslides_med': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - Exposure Index/NOAA_tsunami_aware/landslides_med.shp",
            'rank': 2,
        },
        'landslides_high': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - Exposure Index/NOAA_tsunami_aware/landslides_high.shp",
            'rank': 3,
        },
        'landslides_vhigh': {
            'path': "D:/NFWF_PhaseIII/Guam/Data - Exposure Index/NOAA_tsunami_aware/landslides_vhigh.shp",
            'rank': 4,
        },
    }
    intermediate_dir = "D:/NFWF_PhaseIII/Guam/Data - Exposure Index/threat_index_revisions"
    boundary_path = "D:/Packages/GU_Wildlife_Inputs_v1_01252021_c4e628/commondata/boundaries/GU_30m_Boundary.shp"
    existing_index_path = "D:/Packages/GU_Threat_Inputs_v1_01252021_71a7fb/commondata/raster_data10/GU_Threat_Index_v1.tif"
    intermediate_index_path = "D:/Packages/GU_Threat_Inputs_v1_01252021_71a7fb/commondata/raster_data10/GU_Threat_Index_v2_intermediate.tif"
    revise_index_with_additions(
        new_input_dict, intermediate_dir, boundary_path, existing_index_path,
        intermediate_index_path)

    # add wave exposure, reclassified in Arc to 1-5
    wave_exp_path = "D:/NFWF_PhaseIII/Guam/Data - digital atlas of Southern Guam/coast-wave-exposure/wave_exposure_5class_proj2.tif"
    raster_info = pygeoprocessing.get_raster_info(existing_index_path)
    input_datatype = raster_info['datatype']
    input_nodata = raster_info['nodata'][0]
    pixel_size = raster_info['pixel_size']
    with tempfile.NamedTemporaryFile(
            prefix='revised_index_not_clipped', suffix='.tif',
            delete=False) as unclipped_raster_file:
        unclipped_raster_path = unclipped_raster_file.name
        sum_path_list = [intermediate_index_path, wave_exp_path]
    
    sum_path_list = [intermediate_index_path, wave_exp_path]
    raster_list_sum(
        sum_path_list, input_nodata, unclipped_raster_path, input_nodata,
        input_datatype, nodata_remove=True)

    # clip their sum to the region
    revised_index_path = "D:/Packages/GU_Threat_Inputs_v1_01252021_71a7fb/commondata/raster_data10/GU_Threat_Index_ALL_v2.tif"
    pygeoprocessing.align_and_resize_raster_stack(
        [unclipped_raster_path], [revised_index_path], ['near'],
        pixel_size, 'intersection', base_vector_path_list=[boundary_path],
        raster_align_index=0)
    
    # clean up
    os.remove(unclipped_raster_path)


if __name__ == "__main__":
    # AS_threat_index_v2()
    # AS_terrestrial_index_revisions()
    GU_marine_index_revisions()
    GU_terrestrial_index_revisions()
    GU_threat_index_revisions()