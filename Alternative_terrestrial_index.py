"""Alternative terrestrial index for NFWF."""
import os
import tempfile
import shutil

import numpy

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing


# study area boundary
# Puerto Rico
# _BOUNDARY_AOI = "D:/Packages/Alternative_Terrestrial_Index/Clipped_Puerto_Rico/PR_30mDepth_Wildlife_WGS84.shp"
# Great Lakes
# _BOUNDARY_AOI = "D:/Packages/Alternative_Terrestrial_Index/Great_Lakes/GL_PrelimBndy.shp"

# nodata value for normalized rasters of Int16 type
_TARGET_NODATA = -1


def raster_list_sum(
        raster_list, input_nodata, target_path, target_nodata,
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
            target_path, gdal.GDT_Int16, target_nodata)

    else:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in raster_list], raster_sum_op,
            target_path, gdal.GDT_Int16, target_nodata)


def normalize(value_raster_path, target_path, aoi_path, weight):
    """Calculate a normalized raster according to values in boundary.

    Args:
        value_raster_path (string): path to raster containing values that
            should be normalized
        target_path (string): path to location where normalized raster should
            be created

    Side effects:
        creates or modifies a raster at the location of ``target_path``

    Returns:
        None

    """
    def normalize_op(raster, min_value, max_value, weight_value):
        """Normalize values inside a raster.

        Calculate normalized values that lie between 0 and 1 according to
        (value - min) / (max - min). Then, multiply these values by 100 and
        return an array of normalized integer values lying between 0 and 100.
        Then, multiply these by `weight` in case the normalized version should
        count for more than 1 relative to other normalized inputs.

        Args:
            raster (numpy.ndarray): values to be normalized
            min_value (float or int): minimum value in the raster
            max_value (float or int): maximum value in the raster
            weight_value (float or int): weight to apply to normalized values

        Returns:
            array of normalized integer values multiplied by weight

        """
        valid_mask = (~numpy.isclose(raster, raster_nodata))
        float_ar = numpy.empty(raster.shape, dtype=numpy.float32)
        float_ar[:] = _TARGET_NODATA
        float_ar[valid_mask] = (
            (raster[valid_mask] - min_value) / (max_value - min_value))
        result = numpy.empty(raster.shape, dtype=numpy.int16)
        result[:] = _TARGET_NODATA
        result[valid_mask] = float_ar[valid_mask] * 100 * weight_value
        return result
    # calculate zonal max and min with zonal statistics inside boundary aoi
    zonal_stats = pygeoprocessing.zonal_statistics(
        (value_raster_path, 1), aoi_path)
    zonal_min = zonal_stats[0]['min']
    zonal_max = zonal_stats[0]['max']
    
    # calculate normalized raster according to zonal max and min
    raster_nodata = pygeoprocessing.get_raster_info(
        value_raster_path)['nodata'][0]
    pygeoprocessing.raster_calculator(
            [(value_raster_path, 1), (zonal_min, 'raw'), (zonal_max, 'raw'),
            (weight, 'raw')],
        normalize_op, target_path, gdal.GDT_Int16, _TARGET_NODATA)


def alternative_index_workflow(
        workspace_dir, raster_input_dict, aoi_path, index_path,
        polygon_input_list=None):
    """Compute the alternative index from raw inputs.

    All inputs, including AOI, must be share coordinate reference system and
    must have roughly equivalent extents. Recommend that inputs are clipped and
    projected in Arc prior to running this script.

    Args:
        workspace_dir (string): path to workspace where intermediate results
            should be created/stored
        raster_input_dict (dict): a nested python dictionary containing info
            about raster-based inputs that should be combined. The keys in the
            index should be the labels for each input; values in the dictionary
            should be dictionaries containing the keys 'path' (path to the
            raster input) and 'weight' (weighting value that is applied to the
            normalized values in this input relative to others). EACH INDEX IS
            INTERPRETED AS HIGH VALUE = GOOD.        
        aoi_path (string): path to boundary of the study area
        index_path (string): path to location where the index should be saved
        polygon_input_list (list): list of paths to polygon inputs that should
            be included. Each of these is assigned a weight of 1.

    Side effects:
        creates or modifies a raster at the location ``index_path``

    Returns:
        None
    
    """
    # ensure that each new input shares spatial reference
    vector_info = pygeoprocessing.get_vector_info(aoi_path)
    destination_proj = osr.SpatialReference()
    destination_proj.ImportFromWkt(vector_info['projection_wkt'])
    problem_list = []
    for new_input in raster_input_dict:
        new_proj = osr.SpatialReference()
        new_proj.ImportFromWkt(
            pygeoprocessing.get_raster_info(
                raster_input_dict[new_input]['path'])['projection_wkt'])
        if (new_proj.IsSame(destination_proj) == 0):
            problem_list.append(new_input)
    if problem_list:
        raise ValueError(
            "Project these to match the AOI: {}".format(problem_list))
        
    intermediate_dir = os.path.join(workspace_dir, 'intermediate')
    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)

    normalized_dir = os.path.join(intermediate_dir, 'normalized')
    if not os.path.exists(normalized_dir):
        os.makedirs(normalized_dir)

    aligned_dir = os.path.join(intermediate_dir, 'aligned')
    if not os.path.exists(aligned_dir):
        os.makedirs(aligned_dir)

    # normalize all raster-based inputs within AOI
    base_raster_path_list = []
    aligned_raster_path_list = []
    for new_input in raster_input_dict:
        value_raster_path = raster_input_dict[new_input]['path']
        try:
            weight = raster_input_dict[new_input]['weight']
        except KeyError:
            weight = 1
        bn = os.path.basename(value_raster_path)
        normalized_path = os.path.join(normalized_dir, bn)
        aligned_path = os.path.join(aligned_dir, bn)
        base_raster_path_list.append(normalized_path)
        aligned_raster_path_list.append(aligned_path)
        if not os.path.exists(normalized_path):
            with tempfile.NamedTemporaryFile(
                    prefix='mask_raster', delete=False, suffix='.tif',
                    dir=normalized_dir) as clipped_raster_file:
                clipped_raster_path = clipped_raster_file.name
            pygeoprocessing.mask_raster(
                (value_raster_path, 1), aoi_path, clipped_raster_path)
            normalize(clipped_raster_path, normalized_path, aoi_path, weight)
            os.remove(clipped_raster_path)
    
    # align and resample normalized rasters, using minimum pixel size of inputs
    pixel_size_list = []
    for new_input in raster_input_dict:
        value_raster_path = raster_input_dict[new_input]['path']
        raster_info = pygeoprocessing.get_raster_info(value_raster_path)
        pixel_size_list.append(raster_info['pixel_size'])
    target_pixel_size = min(pixel_size_list)
    min_pixel_index = pixel_size_list.index(min(pixel_size_list))
    
    if not all([os.path.exists(f) for f in aligned_raster_path_list]):
        pygeoprocessing.align_and_resize_raster_stack(
            base_raster_path_list, aligned_raster_path_list,
            ['near'] * len(base_raster_path_list), target_pixel_size,
            'intersection', raster_align_index=min_pixel_index)
    
    # rasterize polygon inputs
    template_raster_path = aligned_raster_path_list[0]
    if polygon_input_list:
        for vec_path in polygon_input_list:
            target_raster_path = os.path.join(
                aligned_dir, '{}.tif'.format(os.path.basename(vec_path)[:-4]))
            aligned_raster_path_list.append(target_raster_path)
            if not os.path.exists(target_raster_path):
                pygeoprocessing.new_raster_from_base(
                    template_raster_path, target_raster_path, gdal.GDT_Int16,
                    [_TARGET_NODATA], fill_value_list=[_TARGET_NODATA])
                pygeoprocessing.rasterize(
                    vec_path, target_raster_path, burn_values=[100])
    
    # add together
    raster_list_sum(
        aligned_raster_path_list, _TARGET_NODATA, index_path, _TARGET_NODATA,
        nodata_remove=True)


def puerto_rico_demonstration():
    """Calculate alternative index for Puerto Rico."""
    workspace_dir = 'D:/Packages/Alternative_Terrestrial_Index/Clipped_Puerto_Rico'
    raster_input_list = [
        os.path.join(workspace_dir, p) for p in [
            'mammals_all_threats.tif', 'Amphibians_all_threats.tif',
            'Biodiversity_Significance.tif', 'Forest_Landscape_Integrity.tif',
            'Landscape_Biodiversity_Intactness_Index.tif']]
    
    # all polygons inside "PrimaryHabitat_Inputs" for Puerto Rico for all taxa
    polygon_input_list = [
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Terrestrial_Mammals/Intermediate_Data/PrimaryHabitat_Inputs/PR_DRNA_CritWildlife_Mammals.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Terrestrial_Mammals/Intermediate_Data/PrimaryHabitat_Inputs/CamuyCave_Bndy.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Reptiles/Intermediate_Data/PrimaryHabitat_Inputs/Reptiles_PRDRNA_CriticalWildlifeAreas.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Reptiles/Intermediate_Data/PrimaryHabitat_Inputs/Reptiles_NOAA_EnvSensitivityIndex.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Reptiles/Intermediate_Data/PrimaryHabitat_Inputs/Reptiles_FederalDesignations.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Amphibians/Intermediate_Data/PrimaryHabitat_Inputs/Amphibian_PRDRNA_CriticalWildlifeAreas.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Amphibians/Intermediate_Data/PrimaryHabitat_Inputs/Amphibian_NOAA_EnvSensitivityIndex.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Amphibians/Intermediate_Data/PrimaryHabitat_Inputs/Amphibian_FederalDesignations.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Birds/Intermediate_Data/PrimaryHabitat_Inputs/Bird_NOAA_EnvSensitivityIndex.shp',
        'D:/NFWF_PhaseII/PRVI/FINAL/_Puerto_Rico/Wildlife_Index/Terrestrial/Birds/Intermediate_Data/PrimaryHabitat_Inputs/Bird_FederalDesignations.shp']
    index_path = os.path.join(
        workspace_dir, 'alternative_index_with_polygons_5.10.21.tif')
    alternative_index_workflow(
        workspace_dir, raster_input_list, polygon_input_list, index_path)


def great_lakes_demonstration():
    """Calculate alternative index for Great Lakes."""
    workspace_dir = 'D:/Packages/Alternative_Terrestrial_Index/Great_Lakes'
    raster_input_list = [
        "D:/Packages/Alternative_Terrestrial_Index/Great_Lakes/Local_Biodiversity_Intactness_Index.tif",
        "D:/Packages/Alternative_Terrestrial_Index/Great_Lakes/Forest_Landscape_Integrity.tif",
        "D:/Packages/Alternative_Terrestrial_Index/Great_Lakes/GAP_habitat_richness_sum.tif"]

    polygon_input_list = [
        "D:/Packages/Alternative_Terrestrial_Index/Great_Lakes/IBAs_20150702.shp",
        "D:/Packages/Alternative_Terrestrial_Index/Great_Lakes/CRITHAB_POLY.shp",
        "D:/Packages/Alternative_Terrestrial_Index/Great_Lakes/PADUS_GAP_1_2.shp"]

    index_path = os.path.join(
        workspace_dir, 'alternative_index_with_polygons_5-14-21.tif')
    alternative_index_workflow(
        workspace_dir, raster_input_list, polygon_input_list, index_path)


def alaska_demonstration():
    """Add FLI and LBII to existing index."""
    workspace_dir = 'D:/Packages/Alternative_Terrestrial_Index/Alaska/index_v1_rank1'
    raster_input_dict = {
        'index_v1': {
            'path': "D:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data23/AK_Terrestrial_Index_All_v1.tif",
            # 'weight': 3,
        },
        'lbii': {
            'path': "D:/Packages/Alternative_Terrestrial_Index/Alaska/Local_Biodiversity_Intactness_proj.tif",
        },
        'fli': {
            'path': "D:/Packages/Alternative_Terrestrial_Index/Alaska/Forest_Landscape_Integrity_proj.tif",
        }
    }

    aoi_path = "D:/Packages/Alternative_Terrestrial_Index/Alaska/AK_20mDepth_Boundary_v1_diss.shp"
    index_path = os.path.join(
        workspace_dir, 'alternative_index_v1_rank1.tif')
    alternative_index_workflow(
        workspace_dir, raster_input_dict, aoi_path, index_path)


if __name__ == "__main__":
    # puerto_rico_demonstration()
    # great_lakes_demonstration()
    alaska_demonstration()
