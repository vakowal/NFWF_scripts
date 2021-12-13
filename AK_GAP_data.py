"""Process GAP Analysis rasters for Alaska."""
import os
import zipfile

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing

# folder containing zipped folders downloaded from GAP Analysis Program
_GAP_FOLDER = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp"

# study area boundary
_BOUNDARY_PATH = "D:/Packages/AK_Wildlife_012721_4c7e75/commondata/boundaries/AK_20mDepth_Boundary_v1.shp"

_TARGET_NODATA = 255

_TARGET_DATATYPE = gdal.GDT_Byte

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

        
def gap_data_workflow():
    # unzip all raw folders
    folder_list = [
        f for f in os.listdir(_GAP_FOLDER) if f.endswith('.zip')]

    # TODO within species group (need a table linking folder name to spp group)
    # dictionary: key is group, value is list of
    # TODO destination folder
    spp_group_dict = {
        'test': [
            'AmericanBeaver_AnnualDistribution',
            'AmericanBlackBear_AnnualDistribution',
            'AmericanMarten_AnnualDistribution'],
        }
    destination = "C:/Users/Ginger/Desktop/test_extract"

    for sp_key in spp_group_dict:
        for sp in spp_group_dict[sp_key]:
            zip_path = os.path.join(_GAP_FOLDER, '{}.zip'.format(sp))
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(destination)

        raster_list = [
            os.path.join(destination, '{}.img'.format(sp)) for sp in
            spp_group_dict[sp_key]]
        input_nodata = pygeoprocessing.get_raster_info(
            raster_list[0])['nodata'][0]
        for raster_path in raster_list:
            test_nodata = pygeoprocessing.get_raster_info(
                raster_path)['nodata'][0]
            if test_nodata != input_nodata:
                raise ValueError("Input rasters do not share nodata value")

        # TODO align rasters in raster_list into temp directory
        # TODO does this work?
        target_path = os.path.join(destination, '{}.tif'.format(sp_key))
        raster_list_sum(
            raster_list, input_nodata, target_path, _TARGET_NODATA,
            _TARGET_DATATYPE, nodata_remove=True)

        # TODO clip to region
        # TODO delete temporary directory


if __name__ == "__main__":
    gap_data_workflow()