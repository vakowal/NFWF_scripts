# Revise threat index for Alaska by using elevation mask
# to decrease floodprone areas
import os
import tempfile
import shutil

import numpy

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing

# all outputs should align with this template raster
_TEMPLATE_PATH = "D:/Current/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data11/AK_Threat_Index_10class_v1.tif"

# items to be added together should have this nodata value
_TARGET_NODATA = 255


def align_inputs(align_dir):
    """Align all inputs to the threat index raster.

    Align all inputs to the threat index raster, storing aligned rasters in
    `align_dir`. Reclassify nodata areas to _TARGET_NODATA.

    Return:
        a dictionary where keys identify inputs and values are paths to
            aligned rasters for each input

    """
    template_raster_info = pygeoprocessing.get_raster_info(_TEMPLATE_PATH)
    base_path_id_map = {
        'low_lying': "D:/Current/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data9/AK_Low_Lying_Areas_v1.tif",
        'erodibility': "D:/Current/Packages/AK_Threat_Inputs_012721_9522b9/AK_Erosion_STA_Add_reclass_revised.tif",
        'permafrost': "D:/Current/Packages/AK_Threat_Inputs_012721_9522b9/AK_Permafrost_STA_Add_revised.tif",
        'tsunami': "D:/Current/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data7/AK_Tsunami_v1.tif",
        'floodprone_no_sta': "D:/Current/Packages/AK_Threat_Inputs_012721_9522b9/commondata/raster_data3/AK_Floodprone_Inputs_v1.tif",
        # revised to include only communities not in group 3 according to STA
        'sta': "D:/Current/Packages/AK_Threat_Inputs_012721_9522b9/commondata/floodprone/AK_STA_Flooding_revised.tif",  
        'dem': "D:/Current/Alaska/Data/ifsar_dem_resample_30m/ifsar_30m_proj.tif", # ifsar DEM resampled to 30 m
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
        print("\nAligning inputs .....")
        pygeoprocessing.align_and_resize_raster_stack(
            base_input_path_list, aligned_path_list,
            ['near'] * len(aligned_path_list),
            template_raster_info['pixel_size'],
            bounding_box_mode=template_raster_info['bounding_box'],
            raster_align_index=0)
        for key in [
                'floodprone_no_sta', 'low_lying', 'erodibility', 'permafrost',
                'tsunami']:
            print("\nReclassifying nodata .....")
            reclassify_nodata(aligned_inputs[key], _TARGET_NODATA)
    return aligned_inputs


def reclassify_nodata(target_path, new_nodata_value):
    """Reclassify the nodata value of a raster to a new value.

    Convert all areas of nodata in the target raster to the new nodata
    value, which must be an integer.

    Parameters:
        target_path (string): path to target raster
        new_nodata_value (integer): new value to set as nodata

    Side effects:
        modifies the raster indicated by `target_path`

    Returns:
        None

    """
    def reclassify_op(target_raster):
        reclassified_raster = numpy.copy(target_raster)
        reclassify_mask = (target_raster == previous_nodata_value)
        reclassified_raster[reclassify_mask] = new_nodata_value
        return reclassified_raster

    fd, temp_path = tempfile.mkstemp()
    shutil.copyfile(target_path, temp_path)
    previous_nodata_value = pygeoprocessing.get_raster_info(
        target_path)['nodata'][0]

    pygeoprocessing.raster_calculator(
        [(temp_path, 1)], reclassify_op, target_path, gdal.GDT_Byte,
        new_nodata_value)

    # clean up
    os.close(fd)
    os.remove(temp_path)


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
            target_path, gdal.GDT_Byte, target_nodata)

    else:
        pygeoprocessing.raster_calculator(
            [(path, 1) for path in raster_list], raster_sum_op,
            target_path, gdal.GDT_Byte, target_nodata)


def revise_floodprone_input(aligned_inputs, elevation_cutoff, target_path):
    """Set areas in floodprone input to 0 where elevation is > cutoff value.

    According to the DEM, set areas in floodprone input to 0 where elevation is
    greater than `elevation_cutoff`. Then mosaic floodprone communities
    identified by statewide threat assessment (STA) to make sure that these are
    included in floodprone input no matter the elevation values.

    Parameters:
        aligned_inputs (dict): dictionary where values are paths to aligned
            rasters, including dem, floodprone areas original, and sta
            floodprone communities
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
            (~numpy.isclose(floodprone_ar, _TARGET_NODATA)) &
            (~numpy.isclose(sta_ar, sta_nodata)))
        zero_mask = (
            (dem_ar > elevation_cutoff) &
            (sta_ar == 0) &
            valid_mask)
        result = numpy.copy(floodprone_ar)
        result[zero_mask] = 0
        return result


    dem_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['dem'])['nodata'][0]
    sta_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['sta'])['nodata'][0]

    # calculate modified floodprone input
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in
            [aligned_inputs['dem'], aligned_inputs['floodprone_orig'],
            aligned_inputs['sta']]],
        elevation_cutoff_op, target_path, gdal.GDT_Byte, _TARGET_NODATA)


def revise_floodprone_input_ifsar(aligned_inputs, target_path):
    """Set areas in floodprone input to 0 where ifsar DEM is nodata.

    The ifsar DEM contains valid values only where elevation is between 0 and
    20 m. Areas in floodprone input should be set to 0 outside this range,
    where the ifsar DEM is nodata. Then mosaic floodprone communities
    identified by statewide threat assessment (STA) to make sure that these are
    included in floodprone input no matter the elevation values.

    Parameters:
        aligned_inputs (dict): dictionary where values are paths to aligned
            rasters, including dem, floodprone areas without sta, and sta
            floodprone communities
        target_path (string): path to save the revised floodprone input

    Returns:
        None

    """
    def elevation_cutoff_op(dem_ar, floodprone_ar, sta_ar):
        """Restrict floodprone to valid areas in ifsar DEM, and add STA."""
        floodprone_revised = numpy.copy(floodprone_ar)

        # set areas outside valid elevation range to 0
        zero_mask = (numpy.isclose(dem_ar, dem_nodata))
        floodprone_revised[zero_mask] = 0

        # mosaic STA values into index
        sta_mask = (sta_ar == 1)
        floodprone_revised[sta_mask] = floodprone_revised[sta_mask] + 1
        return floodprone_revised

    dem_nodata = pygeoprocessing.get_raster_info(
        aligned_inputs['dem'])['nodata'][0]

    # calculate modified floodprone input
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in
            [aligned_inputs['dem'], aligned_inputs['floodprone_no_sta'],
            aligned_inputs['sta']]],
        elevation_cutoff_op, target_path, gdal.GDT_Byte, _TARGET_NODATA)


def revise_threat_index(aligned_inputs, revised_floodprone_path, target_path):
    """Calculate revised threat index that includes revised floodprone areas.

    Parameters:
        aligned_inputs (dict): dictionary where values are paths to aligned
            rasters, including low-lying areas, erodibility, permafrost,
            tsunami risk
        revised_floodprone_path (string): path to revised floodprone areas
        target_path (string): location where revised threat index should be
            saved

    Returns:
        None

    """
    input_key_list = ['low_lying', 'erodibility', 'permafrost', 'tsunami']
    input_path_list = [
        revised_floodprone_path] + [aligned_inputs[v] for v in input_key_list]

    raster_list_sum(
        input_path_list, _TARGET_NODATA, target_path, _TARGET_NODATA)


def threat_revisions_workflow():
    """Calculate revised floodprone areas and revised threat index."""
    output_dir = "E:/Current/Alaska/Revise_threat_index/AK_threat_revise"
    elevation_cutoff_list = [20]  # [25]  # [300, 500, 700]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    align_dir = os.path.join(output_dir, 'aligned_inputs')  # tempfile.mkdtemp()
    aligned_inputs = align_inputs(align_dir)
    for elevation_cutoff in elevation_cutoff_list:
        revised_floodprone_path = os.path.join(
            output_dir, 'floodprone_revised_{}m.tif'.format(elevation_cutoff))
        if not os.path.isfile(revised_floodprone_path):
            revise_floodprone_input(
                aligned_inputs, elevation_cutoff, revised_floodprone_path)
        revised_threat_path = os.path.join(
            output_dir,
            'Threat_Index_revised_{}m.tif'.format(elevation_cutoff))
        if not os.path.isfile(revised_threat_path):
            revise_threat_index(
                aligned_inputs, revised_floodprone_path, revised_threat_path)


def threat_revisions_ifsar_DEM_workflow():
    """Calculate revised floodprone areas and threat index using ifsar DEM.

    This workflow is a little different because the ifsar DEM only contains
    valid values in the range [0, 20]. The application of the elevation mask
    is simpler: only areas with valid values in the DEM should have a valid
    value in the floodprone input to the threat index.

    Returns:
        None

    """
    output_dir = "D:/Current/Alaska/Revise_threat_index/AK_threat_revise_STA_communities_ifsar_dem"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    align_dir = os.path.join(output_dir, 'aligned_inputs')
    if not os.path.exists(align_dir):
        os.makedirs(align_dir)
    aligned_inputs = align_inputs(align_dir)
    revised_floodprone_path = os.path.join(
        output_dir, 'floodprone_revised_ifsar_20m.tif')
    if not os.path.isfile(revised_floodprone_path):
        print("\nCalculating revised floodprone input .....")
        revise_floodprone_input_ifsar(
            aligned_inputs, revised_floodprone_path)
    revised_threat_path = os.path.join(
        output_dir, 'Threat_Index_revised_ifsar_20m.tif')
    if not os.path.isfile(revised_threat_path):
        print("\nCalculating revised threat index .....")
        revise_threat_index(
            aligned_inputs, revised_floodprone_path, revised_threat_path)


def mosaic_zonal_mean_with_surface():
    """Mosaic community footprints with continuous index.

    Where community footprints have values >0, they should supersede values in
    the continuous index.

    Returns:
        None

    """
    def mosaic_op(footprint_ar, index_ar):
        """Mosaic footprints into continuous index."""
        mosaic_mask = (footprint_ar > 0)
        result = numpy.copy(index_ar)
        result[mosaic_mask] = footprint_ar[mosaic_mask]
        return result

    # outer_dir = "D:/Current/Packages/AK_Threat_Inputs_012721_9522b9"
    outer_dir = "D:/Current/Packages/AK_Asset_Inputs_012721_2fd5ef"

    # align index and community footprints together
    align_dir = tempfile.mkdtemp()
    base_path_id_map = {
        # 'footprints': os.path.join(
        #     outer_dir, "AK_Communities_Subtract3_Threat_mean_int.tif"),
        # 'threat_index': os.path.join(
        #     outer_dir, "AK_Threat_Index_10class_revised.tif"),
        'footprints': os.path.join(
            outer_dir, 'AK_STA_Communities_Exposure_mean_int.tif'),
        'index': os.path.join(
            outer_dir, 'AK_Exposure_Index_10class_revised1.tif'),
    }
    base_input_path_list = [
        base_path_id_map[k] for k in sorted(base_path_id_map.keys())]
    aligned_inputs = dict([(key, os.path.join(
        align_dir, 'aligned_%s' % os.path.basename(path)))
        for key, path in base_path_id_map.items()])
    aligned_path_list = [
        aligned_inputs[k] for k in sorted(aligned_inputs.keys())]
    
    template_raster_info = pygeoprocessing.get_raster_info(
        base_path_id_map['index'])
    pygeoprocessing.align_and_resize_raster_stack(
        base_input_path_list, aligned_path_list,
        ['near'] * len(aligned_path_list),
        template_raster_info['pixel_size'],
        bounding_box_mode=template_raster_info['bounding_box'])

    # mosaic community footprints into index raster
    # target_path = os.path.join(
    #     outer_dir, 'AK_Threat_Index_10class_revised_footprint_mean.tif')
    target_path = os.path.join(
        outer_dir, 'AK_Exposure_Index_10class_revised1_footprint_mean.tif')
    pygeoprocessing.raster_calculator(
        [(path, 1) for path in
            [aligned_inputs['footprints'], aligned_inputs['index']]],
        mosaic_op, target_path, gdal.GDT_Byte,
        template_raster_info['nodata'][0])

    # clean up
    shutil.rmtree(align_dir)


if __name__ == "__main__":
    # threat_revisions_workflow()
    # threat_revisions_ifsar_DEM_workflow()
    mosaic_zonal_mean_with_surface()
