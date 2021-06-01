"""Calculate threat index v2 for American Samoa."""
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


def calc_threat_index_v2():
    # threat index v1, all values (not reclassified)
    threat_index_v1_path = "D:/Packages/AS_Threat_Inputs_v1_01252021_de5873/commondata/raster_data10/AS_Threat_Index_v1.tif"
    # polygon data showing extent of tsunami inundation
    tsunami_shp_path = "D:/NFWF_PhaseIII/American_Samoa/Data - Exposure Index/Threat_Index/Tsunami/AmericanSamoa_PMTInundationZone_slowWalk_110_timeMap_proj.shp"
    # NFWF boundary
    reg_shp_path = "D:/Packages/AS_Threat_Inputs_v1_01252021_de5873/commondata/boundaries/AS_30mDepth_Bndy.shp"
    tsunami_raster_path = "D:/NFWF_PhaseIII/American_Samoa/Data - Exposure Index/Threat_Index/Tsunami/AmericanSamoa_PMTInundationZone_slowWalk_110_timeMap_val3.tif"
    threat_index_v2_path = "D:/Packages/AS_Threat_Inputs_v1_01252021_de5873/commondata/raster_data10/AS_Threat_Index_v2.tif"
    raster_info = pygeoprocessing.get_raster_info(threat_index_v1_path)
    input_datatype = raster_info['datatype']
    input_nodata = raster_info['nodata'][0]
    pixel_size = raster_info['pixel_size']
    tsunami_val = 3
    # rasterize tsunami, using v1 threat index as template
    # pygeoprocessing.new_raster_from_base(
    #     threat_index_v1_path, tsunami_raster_path, input_datatype,
    #     [input_nodata], fill_value_list=[input_nodata])
    # pygeoprocessing.rasterize(
    #     tsunami_shp_path, tsunami_raster_path, burn_values=[tsunami_val])
    # add the two together; nodata reclassified to zero
    with tempfile.NamedTemporaryFile(
            prefix='threat_v2_not_clipped', suffix='.tif',
            delete=False) as unclipped_raster_file:
        unclipped_raster_path = unclipped_raster_file.name
    raster_list_sum(
        [threat_index_v1_path, tsunami_raster_path], input_nodata,
        unclipped_raster_path, input_nodata, input_datatype,
        nodata_remove=True)
    # clip their sum to the region
    pygeoprocessing.align_and_resize_raster_stack(
        [unclipped_raster_path], [threat_index_v2_path], ['near'],
        pixel_size, 'intersection', base_vector_path_list=[reg_shp_path],
        raster_align_index=0)

    # clean up
    os.remove(unclipped_raster_path)


if __name__ == "__main__":
    calc_threat_index_v2()