"""Prepare geotiffs for CREST."""
import os
import tempfile
import shutil

import numpy

from osgeo import gdal

import pygeoprocessing

# rasters for CREST should have this data type
_TARGET_DATATYPE = gdal.GDT_Byte

# rasters for CREST should have this nodata value
_TARGET_NODATA = 255

# boundary for Guam
_BOUNDARY_SHP = "D:/Packages/GU_Wildlife_Inputs_v1_01252021_c4e628/commondata/boundaries/GU_30m_Boundary.shp"

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
        print('processing raster: {}'.format(bn))
        reclassify_nodata(input_path, target_path)


def check_guam_rasters():
    """double check rasters uploaded to CREST."""
    target_raster_dir = "D:/NFWF_PhaseIII/Guam/FOR CREST/for_upload"
    input_path_list = [
        f for f in os.listdir(target_raster_dir) if f.endswith('.tif')]
    for input_bn in input_path_list:
        input_path = os.path.join(target_raster_dir, input_bn)
        raster_info = pygeoprocessing.get_raster_info(input_path)
        pixel_size = raster_info['pixel_size'][0]
        nodata = raster_info['nodata'][0]
        datatype = raster_info['datatype']
        max_val = pygeoprocessing.zonal_statistics(
            (input_path, 1), _BOUNDARY_SHP)[0]['max']
        print("stats for {}:".format(os.path.basename(input_path)))
        print("datatype: {}, pixelsize: {}, nodata: {}, max val: {}".format(
            datatype, pixel_size, nodata, max_val))


def resample_guam_rasters():
    """Align a few bad rasters to identical pixel size."""
    template_raster = "D:/NFWF_PhaseIII/Guam/FOR CREST/for_upload/GU_Threat_Index_10class_v2.tif"
    pixel_size = pygeoprocessing.get_raster_info(template_raster)['pixel_size']
    with tempfile.NamedTemporaryFile(
            prefix='template_aligned', suffix='.tif',
            delete=False) as template_aligned_file:
        template_aligned_path = template_aligned_file.name
    raw_input_list = [
        template_raster,
        "D:/NFWF_PhaseIII/Guam/FOR CREST/raw_rasters_from_mxd/GU_Slope_v1.tif",
        "D:/NFWF_PhaseIII/Guam/FOR CREST/raw_rasters_from_mxd/GU_Terrestrial_index_4class_v3.tif",
        "D:/NFWF_PhaseIII/Guam/FOR CREST/raw_rasters_from_mxd/wave_exposure_5class_proj.tif"]
    aligned_path_list = [
        template_aligned_path,
        "D:/NFWF_PhaseIII/Guam/FOR CREST/for_upload/GU_Slope_v1.tif",
        "D:/NFWF_PhaseIII/Guam/FOR CREST/for_upload/GU_Terrestrial_index_4class_v3.tif",
        "D:/NFWF_PhaseIII/Guam/FOR CREST/for_upload/wave_exposure_5class_proj.tif"]
    pygeoprocessing.align_and_resize_raster_stack(
        raw_input_list, aligned_path_list, ['near'] * len(raw_input_list),
        pixel_size, 'union')
    check_guam_rasters()


if __name__ == "__main__":
    # process_guam_rasters()
    # check_guam_rasters()
    # resample_guam_rasters()
    
