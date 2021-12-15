"""Process GAP Analysis rasters for Alaska."""
import os
import zipfile
import tempfile
import shutil

import numpy

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing

# folder containing zipped folders downloaded from GAP Analysis Program
_GAP_FOLDER = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp"

# study area boundary
_BOUNDARY_PATH = "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/boundaries/AK_20mDepth_Boundary_v1.shp"

_TARGET_NODATA = 255

_TARGET_DATATYPE = gdal.GDT_Byte


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
    previous_raster_info = pygeoprocessing.get_raster_info(target_path)
    previous_nodata_value = previous_raster_info['nodata'][0]
    datatype = previous_raster_info['datatype']

    pygeoprocessing.raster_calculator(
        [(temp_path, 1)], reclassify_op, target_path, datatype,
        new_nodata_value)

    # clean up
    os.close(fd)
    os.remove(temp_path)


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
    """Process GAP data for Alaska subsistence species.

    Unzip folders downloaded from the Alaska GAP analysis project. Within
    species groups, sum up the distribution rasters for individual species.
    Clip the sum of the species distribution rasters inside the species group
    to the study area boundary.

    Returns:
        None

    """
    # species groups
    spp_group_dict = {
        'Fowl': [
            'AmericanWigeon_BreedingDistribution',
            'ArcticTern_BreedingDistribution',
            'BlueWingedTeal_BreedingDistribution',
            'Brant_BreedingDistribution',
            'Bufflehead_BreedingDistribution',
            'Canvasback_BreedingDistribution',
            'CommonGoldeneye_BreedingDistribution',
            'CommonMerganser_BreedingDistribution',
            'CommonSnipe_BreedingDistribution',
            'Gadwall_AnnualDistribution',
            'GreenWingedTeal_BreedingDistribution',
            'HarlequinDuck_BreedingDistribution',
            'LeastAuklet_BreedingDistribution',
            'LongTailedDuck_BreedingDistribution',
            'Mallard_BreedingDistribution',
            'MarbledMurrelet_BreedingDistribution',
            'NorthernPintail_BreedingDistribution',
            'NorthernShoveler_BreedingDistribution',
            'RedBreastedMerganser_BreedingDistribution',
            'RockPtarmigan_AnnualDistribution',
            'SandhillCrane_BreedingDistribution',
            'SharpTailedGrouse_AnnualDistribution',
            'SnowGoose_BreedingDistribution',
            'SpruceGrouse_AnnualDistribution',
            'SurfScoter_BreedingDistribution',
            'TundraSwan_BreedingDistribution',
            'WillowPtarmigan_AnnualDistribution'],
        'Furbearers': [
            'AmericanBeaver_AnnualDistribution',
            'AmericanMarten_AnnualDistribution',
            'AmericanMink_AnnualDistribution',
            'CanadianLynx_AnnualDistribution',
            'Ermine_AnnualDistribution',
            'LeastWeasel_AnnualDistribution',
            'Muskrat_AnnualDistribution',
            'NorthAmericanPorcupine_AnnualDistribution',
            'RedFox_AnnualDistribution',
            'SnowshoeHare_AnnualDistribution',
            'Wolverine_AnnualDistribution',
            'Woodchuck_AnnualDistribution'],
        'Large_Land_Mammals': [
            'AmericanBlackBear_AnnualDistribution',
            'BrownBear_AnnualDistribution',
            'Coyote_AnnualDistribution',
            'DallsSheep_AnnualDistribution',
            'Moose_AnnualDistribution',
            'MountainGoat_AnnualDistribution',
            'Muskox_AnnualDistribution',
            'Wolf_AnnualDistribution'],
        }
    extract_folder = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Unzip"
    destination_folder = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Processed"

    for sp_key in spp_group_dict:
        print("** working on species group {} **".format(sp_key))
        with tempfile.TemporaryDirectory(
                prefix=sp_key, dir=extract_folder) as temp_spp_dir:
            for sp in spp_group_dict[sp_key]:
                zip_path = os.path.join(_GAP_FOLDER, '{}.zip'.format(sp))
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_folder)

            print("testing nodata values")
            raster_list = [
                os.path.join(extract_folder, '{}.img'.format(sp)) for sp in
                spp_group_dict[sp_key]]
            input_nodata = pygeoprocessing.get_raster_info(
                raster_list[0])['nodata'][0]
            pixel_size = pygeoprocessing.get_raster_info(
                raster_list[0])['pixel_size']
            for raster_path in raster_list:
                test_nodata = pygeoprocessing.get_raster_info(
                    raster_path)['nodata'][0]
                if test_nodata != input_nodata:
                    reclassify_nodata(raster_path, input_nodata)

            print("aligning input rasters")
            aligned_path_list = [
                os.path.join(temp_spp_dir, '{}.tif'.format(sp)) for sp in
                spp_group_dict[sp_key]]
            pygeoprocessing.align_and_resize_raster_stack(
                raster_list, aligned_path_list, ['near'] * len(raster_list),
                pixel_size, 'union')

            print("calculating sum")
            intermediate_path = os.path.join(
                temp_spp_dir, '{}.tif'.format(sp_key))
            raster_list_sum(
                aligned_path_list, input_nodata, intermediate_path,
                _TARGET_NODATA, _TARGET_DATATYPE, nodata_remove=True)

            print("clipping to region")
            target_path = os.path.join(
                destination_folder, '{}.tif'.format(sp_key))
            pygeoprocessing.align_and_resize_raster_stack(
                [intermediate_path], [target_path], ['near'], pixel_size,
                'intersection', base_vector_path_list=[_BOUNDARY_PATH],
                raster_align_index=0,
                vector_mask_options={'mask_vector_path': _BOUNDARY_PATH})


def add_elk_to_large_land_mammals():
    """Add the elk raster to the sum of other large land mammals."""
    extract_folder = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Unzip"
    with tempfile.TemporaryDirectory(dir=extract_folder) as temp_spp_dir:
        raster_list = [
            "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Processed/Large_Land_mammals.tif",
            "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Processed/Elk_AHMG.tif"]
        input_nodata = pygeoprocessing.get_raster_info(
            raster_list[0])['nodata'][0]
        pixel_size = pygeoprocessing.get_raster_info(
            raster_list[0])['pixel_size']
        for raster_path in raster_list:
            test_nodata = pygeoprocessing.get_raster_info(
                raster_path)['nodata'][0]
            if test_nodata != input_nodata:
                reclassify_nodata(raster_path, input_nodata)

        aligned_path_list = [
            os.path.join(temp_spp_dir, os.path.basename(f)) for f in
            raster_list]
        pygeoprocessing.align_and_resize_raster_stack(
            raster_list, aligned_path_list, ['near'] * len(raster_list),
            pixel_size, 'union')

        destination_path = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Processed/Large_Land_Mammals_with_Elk.tif"
        raster_list_sum(
            aligned_path_list, input_nodata, destination_path,
            _TARGET_NODATA, _TARGET_DATATYPE, nodata_remove=True)


if __name__ == "__main__":
    # gap_data_workflow()
    add_elk_to_large_land_mammals()