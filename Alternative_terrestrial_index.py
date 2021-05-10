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
_BOUNDARY_AOI = "D:/Packages/Alternative_Terrestrial_Index/Clipped_Puerto_Rico/PR_30mDepth_Wildlife_WGS84.shp"

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


def normalize(value_raster_path, target_path):
	"""Calculate a normalized raster according to values in _BOUNDARY_AOI.

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
	def normalize_op(raster, min_value, max_value):
		"""Normalize values inside a raster.

		Calculate normalized values that lie between 0 and 1 according to
		(value - min) / (max - min). Then, multiply these values by 100 and
		return an array of normalized integer values lying between 0 and 100.

		Args:
			raster (numpy.ndarray): values to be normalized
			min_value (float or int): minimum value in the raster
			max_value (float or int): maximum value in the raster

		Returns:
			array of normalized integer values between 0 and 100

		"""
		valid_mask = (~numpy.isclose(raster, raster_nodata))
		float_ar = numpy.empty(raster.shape, dtype=numpy.float32)
		float_ar[:] = _TARGET_NODATA
		float_ar[valid_mask] = (
			(raster[valid_mask] - min_value) / (max_value - min_value))
		result = numpy.empty(raster.shape, dtype=numpy.int16)
		result[:] = _TARGET_NODATA
		result[valid_mask] = float_ar[valid_mask] * 100
		return result
	# calculate zonal max and min with zonal statistics inside boundary aoi
	zonal_stats = pygeoprocessing.zonal_statistics(
		(value_raster_path, 1), _BOUNDARY_AOI)
	zonal_min = zonal_stats[0]['min']
	zonal_max = zonal_stats[0]['max']
	
	# calculate normalized raster according to zonal max and min
	raster_nodata = pygeoprocessing.get_raster_info(
		value_raster_path)['nodata'][0]
	pygeoprocessing.raster_calculator(
		[(value_raster_path, 1), (zonal_min, 'raw'), (zonal_max, 'raw')],
		normalize_op, target_path, gdal.GDT_Int16, _TARGET_NODATA)


def alternative_index_workflow(
		workspace_dir, input_bn_list, polygon_input_list, index_path):
	"""Compute the alternative index from raw inputs.

	Args:
		workspace_dir (string): path to workspace where raster-based inputs
			can be found
		input_bn_list (list): list of base names for raster-based inputs that
			should be combined. EACH INDEX IS INTERPRETED AS HIGH VALUE = GOOD
		polygon_input_list (list): list of paths to polygon inputs that should
			be included.
		index_path (string): path to location where the index should be saved

	Side effects:
		creates or modifies a raster at the location ``index_path``

	Returns:
		None
	
	"""
	# normalize all raster-based inputs
	intermediate_dir = os.path.join(workspace_dir, 'intermediate')
	if not os.path.exists(intermediate_dir):
		os.makedirs(intermediate_dir)

	normalized_dir = os.path.join(intermediate_dir, 'normalized')
	if not os.path.exists(normalized_dir):
		os.makedirs(normalized_dir)

	aligned_dir = os.path.join(intermediate_dir, 'aligned')
	if not os.path.exists(aligned_dir):
		os.makedirs(aligned_dir)

	base_raster_path_list = []
	aligned_raster_path_list = []
	for bn in input_bn_list:
		value_raster_path = os.path.join(workspace_dir, bn)
		normalized_path = os.path.join(normalized_dir, bn)
		aligned_path = os.path.join(aligned_dir, bn)
		# normalize(value_raster_path, normalized_path)
		base_raster_path_list.append(normalized_path)
		aligned_raster_path_list.append(aligned_path)
	
	# align and resample normalized rasters, using minimum pixel size of inputs
	pixel_size_list = []
	for bn in input_bn_list:
		pixel_size_list.append(
			pygeoprocessing.get_raster_info(
				os.path.join(workspace_dir, bn))['pixel_size'])
	target_pixel_size = min(pixel_size_list)
	min_pixel_index = pixel_size_list.index(min(pixel_size_list))
	
	# pygeoprocessing.align_and_resize_raster_stack(
	# 	base_raster_path_list, aligned_raster_path_list,
	# 	['near'] * len(base_raster_path_list), target_pixel_size, 'union',
	# 	raster_align_index=min_pixel_index)
	
	# rasterize polygon inputs
	template_raster_path = aligned_raster_path_list[0]
	for vec_path in polygon_input_list:
		target_raster_path = os.path.join(
			aligned_dir, os.path.basename(vec_path))
		aligned_raster_path_list.append(target_raster_path)
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
	input_bn_list = [
		'mammals_all_threats.tif', 'Amphibians_all_threats.tif',
		'Biodiversity_Significance.tif', 'Forest_Landscape_Integrity.tif',
		'Landscape_Biodiversity_Intactness_Index.tif']
	
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
		workspace_dir, input_bn_list, polygon_input_list, index_path)


if __name__ == "__main__":
	puerto_rico_demonstration()