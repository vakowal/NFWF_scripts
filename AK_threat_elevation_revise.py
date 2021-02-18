# Revise threat index for Alaska by using elevation mask
# to decrease floodprone areas
import os
import tempfile

import numpy

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing

# all outputs should align with this template raster
_TEMPLATE_PATH = "E:/Current/Alaska/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data11/AK_Threat_Index_10class_v1.tif"


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
            target_path, gdal.GDT_Float32, target_nodata)

    else:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in raster_list], raster_sum_op,
            target_path, gdal.GDT_Float32, target_nodata)


def revise_floodprone_input(elevation_cutoff, target_path):
    """Set areas in floodprone input to 0 where elevation is > cutoff value.

    According to the DEM, set areas in floodprone input to 0 where elevation is
    greater than `elevation_cutoff`. Then mosaic floodprone communities
    identified by statewide threat assessment (STA) to make sure that these are
    included in floodprone input no matter the elevation values.

    Parameters:
        elevation_cutoff (integer): threshold value. Areas with elevation above
            this value are not counted as floodprone, unless they were
            identified as floodprone by the statewide threat assessment
        target_path (string): path to save the revised floodprone input

    Returns:
        None

    """
    def elevation_cutoff_op(dem_ar, floodprone_ar, sta_ar):
        """Set areas > elevation_cutoff to 0 except as identified by STA."""
        valid_mask = (
            (~numpy.isclose(dem_ar, dem_nodata)) &
            (~numpy.isclose(floodprone_ar, floodprone_nodata)) &
            (~numpy.isclose(sta_ar, sta_nodata)))
        zero_mask = (
            (dem_ar > elevation_cutoff) &
            (sta_ar == 0) &
            valid_mask)
        result = numpy.copy(floodprone_ar)
        result[zero_mask] = 0
        return result

    floodprone_orig_path = "E:/Current/Alaska/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data3/AK_Floodprone_Areas_shift_v1.tif"
    dem_path = "E:/Current/Alaska/Data/elevation_30m_resample_clip/ak_elevation_30m_resample_clip.tif"
    sta_floodprone_path = "E:/Current/Alaska/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data3/AK_Community_STA_Flood_shift.tif"

    # align all with template raster
    template_raster_info = pygeoprocessing.get_raster_info(
        _TEMPLATE_PATH)
    align_dir = "C:/Users/Ginger/Desktop/AK_threat_revise/aligned_inputs"  # tempfile.mkdtemp()
    base_path_id_map = {
        'dem': dem_path,
        'floodprone_orig': floodprone_orig_path,
        'sta': sta_floodprone_path,
    }

    base_input_path_list = [
        base_path_id_map[k] for k in sorted(base_path_id_map.keys())]
    base_input_path_list.insert(0, _TEMPLATE_PATH)
    aligned_inputs = dict([(key, os.path.join(
        align_dir, 'aligned_%s' % os.path.basename(path)))
        for key, path in base_path_id_map.items()])
    aligned_path_list = [
        aligned_inputs[k] for k in sorted(aligned_inputs.keys())]
    aligned_path_list.insert(0, os.path.join(align_dir, 'template.tif'))

    # ensure all inputs match projection of template raster
    for input_path in base_input_path_list:
        input_wkt = pygeoprocessing.get_raster_info(
            input_path)['projection_wkt']
        if input_wkt != template_raster_info['projection_wkt']:
            raise ValueError("Inputs must share projection")

    if not all([os.path.isfile(p) for p in aligned_path_list]):
        pygeoprocessing.align_and_resize_raster_stack(
            base_input_path_list, aligned_path_list,
            ['near'] * len(aligned_path_list),
            template_raster_info['pixel_size'],
            bounding_box_mode=template_raster_info['bounding_box'],
            raster_align_index=0)

    dem_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['dem'])['nodata'][0]
    floodprone_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['floodprone_orig'])['nodata'][0]
    sta_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['sta'])['nodata'][0]

    # calculate modified floodprone input
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in
            [aligned_inputs['dem'], aligned_inputs['floodprone_orig'],
            aligned_inputs['sta']]],
        elevation_cutoff_op, target_path, gdal.GDT_Byte, floodprone_nodata)


def revise_threat_index(revised_floodprone_path, target_path):
    """Calculate revised threat index that includes revised floodprone areas"""
    threat_input_dict = {
        'low_lying': "E:/Current/Alaska/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data9/AK_Low_Lying_Areas_v1.tif",
        'erodibility': "E:/Current/Alaska/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data/AK_Soil_Erodibility_shift_v1.tif",
        'permafrost': "E:/Current/Alaska/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data2/AK_Permafrost_STA_Add.tif",
        'tsunami': "E:/Current/Alaska/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data7/AK_Tsunami_v1.tif",
    }
    input_path_list = [
        revised_floodprone_path] + [v for v in threat_input_dict.values()]

    target_nodata = pygeoprocessing.get_raster_info(
        _TEMPLATE_PATH)['nodata'][0]
    raster_list_sum(
        input_path_list, target_nodata, target_path, target_nodata)


if __name__ == "__main__":
    output_dir = "C:/Users/Ginger/Desktop/threat_index_revisions"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for elevation_cutoff in [300, 500, 700]:
        revised_floodprone_path = os.path.join(
            output_dir, 'floodprone_revised_{}m.tif'.format(elevation_cutoff))
        if not os.path.isfile(revised_floodprone_path):
            revise_floodprone_input(elevation_cutoff, revised_floodprone_path)
        revised_threat_path = os.path.join(
            output_dir,
            'Threat_Index_revised_{}m.tif'.format(elevation_cutoff))
        if not os.path.isfile(revised_threat_path):
            revise_threat_index(revised_floodprone_path, revised_threat_path)
