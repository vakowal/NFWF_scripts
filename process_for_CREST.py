"""Prepare geotiffs for CREST."""
import os

import numpy

from osgeo import gdal

import pygeoprocessing

# rasters for CREST should have this data type
_TARGET_DATATYPE = gdal.GDT_Byte

# rasters for CREST should have this nodata value
_TARGET_NODATA = 255


def reclassify_nodata(input_path, target_path):
    """Reclassify the nodata value of a raster to a default value.

    Convert all areas of nodata in the input raster to the default nodata
    value and convert the data type of the input raster to the default data
    type.
    
    Parameters:
        input_path (string): path to input raster that should be reclassified
        target_path (string): path to target raster with reclassified
    Side effects:
        creates a raster at the location indicated by `target_path`
    
    Returns:
        None
    """
    def reclassify_op(target_raster):
        reclassified_raster = numpy.copy(target_raster)
        reclassify_mask = numpy.isclose(target_raster, previous_nodata_value)
        reclassified_raster[reclassify_mask] = _TARGET_NODATA
        return reclassified_raster

    previous_nodata_value = pygeoprocessing.get_raster_info(
        input_path)['nodata'][0]
    pygeoprocessing.raster_calculator(
        [(input_path, 1)], reclassify_op, target_path, _TARGET_DATATYPE,
        _TARGET_NODATA)


def process_guam_rasters():
    """Process rasters for Guam to be uploaded to CREST.

    Requirements:
        - target data type
        - target nodata value
        - lzw compression

    """
    input_raster_dir = "D:/NFWF_PhaseIII/Guam/FOR CREST/raw_rasters_from_mxd"
    target_raster_dir = "D:/NFWF_PhaseIII/Guam/FOR CREST/for_upload"
    if not os.path.exists(target_raster_dir):
        os.makedirs(target_raster_dir)
    input_bn_list = [
        os.path.basename(f) for f in os.listdir(input_raster_dir) if
        f.endswith('.tif')]
    for bn in input_bn_list:
        input_path = os.path.join(input_raster_dir, bn)
        target_path = os.path.join(target_raster_dir, bn)
        reclassify_nodata(input_path, target_path)


if __name__ == "__main__":
    process_guam_rasters()
