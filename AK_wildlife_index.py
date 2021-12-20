"""Revise Fish and Wildlife Index for Alaska."""
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
        raster_list, target_path, target_datatype, nodata_remove=False):
    """Calculate the sum per pixel across rasters in a list.

    Sum the rasters in `raster_list` element-wise, allowing nodata values
    in the rasters to propagate to the result or treating nodata as zero. If
    nodata is treated as zero, areas where all inputs are nodata will be nodata
    in the output.
    This function also checks for agreement among the inputs in terms of
    spatial reference, pixel size, and nodata values; it aligns all of the
    inputs; and it clips the sum of the inputs to the regional boundary.

    Parameters:
        raster_list (list): list of paths to rasters to sum
        target_path (string): path to location to store the result
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
            numpy.isclose(numpy.array(raster_list), _TARGET_NODATA), axis=0)
        for r in raster_list:
            numpy.place(r, numpy.isclose(r, _TARGET_NODATA), [0])
        sum_of_rasters = numpy.sum(raster_list, axis=0)
        sum_of_rasters[invalid_mask] = _TARGET_NODATA
        return sum_of_rasters

    def raster_sum_op_nodata_remove(*raster_list):
        """Add the rasters in raster_list, treating nodata as zero."""
        invalid_mask = numpy.all(
            numpy.isclose(numpy.array(raster_list), _TARGET_NODATA), axis=0)
        for r in raster_list:
            numpy.place(r, numpy.isclose(r, _TARGET_NODATA), [0])
        sum_of_rasters = numpy.sum(raster_list, axis=0)
        sum_of_rasters[invalid_mask] = _TARGET_NODATA
        return sum_of_rasters

    # check that all inputs have same spatial reference, pixel size, & nodata
    raster_info = pygeoprocessing.get_raster_info(raster_list[0])
    destination_proj = osr.SpatialReference()
    destination_proj.ImportFromWkt(raster_info['projection_wkt'])
    srs_problems = []
    pixel_size_problems = []
    for other_path in raster_list:
        other_info = pygeoprocessing.get_raster_info(other_path)
        new_proj = osr.SpatialReference()
        new_proj.ImportFromWkt(other_info['projection_wkt'])
        if (new_proj.IsSame(destination_proj) == 0):
            srs_problems.append(other_path)
        if other_info['pixel_size'] != raster_info['pixel_size']:
            pixel_size_problems.append(other_path)
    if srs_problems:
        srs_problems.append(raster_list[0])
        raise ValueError(
            "Project these to match: {}".format(srs_problems))
    if pixel_size_problems:
        pixel_size_problems.append(raster_list[0])
        # raise ValueError(
        #     "These do not have same pixel size: {}".format(
        #         pixel_size_problems))  # TODO put back in ?
    for other_path in raster_list:
        other_nodata = pygeoprocessing.get_raster_info(other_path)['nodata'][0]
        if other_nodata != _TARGET_NODATA:
            print("reclassifying nodata")
            reclassify_nodata(other_path, _TARGET_NODATA)

    # align all the rasters to be summed together
    print("aligning {} rasters".format(len(raster_list)))
    with tempfile.TemporaryDirectory() as align_dir:
        aligned_path_list = [
            os.path.join(align_dir, os.path.basename(p)) for p in raster_list]
        pygeoprocessing.align_and_resize_raster_stack(
            raster_list, aligned_path_list, ['near'] * len(raster_list),
            raster_info['pixel_size'], 'union')

        print("calculating raster sum")
        unclip_path = os.path.join(align_dir, 'sum_unclip.tif')
        if nodata_remove:
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in aligned_path_list],
                raster_sum_op_nodata_remove, unclip_path, target_datatype,
                _TARGET_NODATA)
        else:
            pygeoprocessing.raster_calculator(
                [(path, 1) for path in aligned_path_list], raster_sum_op,
                unclip_path, target_datatype, _TARGET_NODATA)

        # clip the sum to the regional boundary
        print("clipping to regional boundary")
        pygeoprocessing.mask_raster(
            (unclip_path, 1), _BOUNDARY_PATH, target_path)
    

def species_groups_sum():
    """Add up species distribution rasters within species groups"""
    gap_folder = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Unzip"
    destination_folder = "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2"
    inputs_dir = {
        # 'crabs': [
        #     "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data/Red_King_Crab_extract.tif",
        #     "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data/Dungeness_Crab.tif",
        #     "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data18/Crabs_EFH.tif"],
        # 'felids_canids': [
        #     "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data48/Felids_Canines_Distributions_Add.tif",
        #     os.path.join(gap_folder, "Wolf_AnnualDistribution.img"),
        #     os.path.join(gap_folder, "RedFox_AnnualDistribution.img"), 
        #     os.path.join(gap_folder, "Coyote_AnnualDistribution.img"),
        #     os.path.join(gap_folder, "CanadianLynx_AnnualDistribution.img")],
        # 'furbearers': [
        #     # "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/Regional_boundary.tif",
        #     os.path.join(gap_folder, 'AmericanBeaver_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'Ermine_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'AmericanMarten_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'Woodchuck_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'AmericanMink_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'Muskrat_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'LeastWeasel_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'NorthAmericanPorcupine_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'Wolverine_AnnualDistribution.img')],
        # 'land_mammals': [
        #     "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data47/Caribou.tif",
        #     os.path.join(gap_folder, 'DallsSheep_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'Muskox_AnnualDistribution.img'),
        #     "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/Elk_AHMG.tif",
        #     os.path.join(gap_folder, 'MountainGoat_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'Moose_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'AmericanBlackBear_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'BrownBear_AnnualDistribution.img')],
        # 'land_birds': [
        #     "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data22/Land_Birds_Distributions_Add_v1.tif",
        #     os.path.join(gap_folder, "RockPtarmigan_AnnualDistribution.img"),
        #     os.path.join(gap_folder, "WillowPtarmigan_AnnualDistribution.img"),
        #     os.path.join(gap_folder, "SpruceGrouse_AnnualDistribution.img"),
        #     os.path.join(gap_folder, "SharpTailedGrouse_AnnualDistribution.img"),
        #     os.path.join(gap_folder, "RuffedGrouse_AnnualDistribution.img")], 
        # 'sea_birds': [
        #     "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data25/Sea_Birds_Distributions_Add_v1.tif",
        #     os.path.join(gap_folder, 'AmericanWigeon_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'ArcticTern_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'BlueWingedTeal_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'Brant_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'Bufflehead_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'Canvasback_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'CommonGoldeneye_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'CommonMerganser_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'CommonSnipe_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'Gadwall_AnnualDistribution.img'),
        #     os.path.join(gap_folder, 'GreenWingedTeal_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'HarlequinDuck_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'LeastAuklet_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'LongTailedDuck_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'Mallard_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'MarbledMurrelet_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'NorthernPintail_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'NorthernShoveler_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'RedBreastedMerganser_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'SandhillCrane_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'SnowGoose_BreedingDistribution.img'),
        #     os.path.join(gap_folder, 'SurfScoter_BreedingDistribution.img')],
        'marine_inverts': [
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/Regional_boundary.tif",
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/blue_mussels_kbec_inverts.tif",
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/cockles.tif",
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/Razor_clam.tif"],
    }
    for spp_group in inputs_dir:
        target_path = os.path.join(
            destination_folder, '{}_sum.tif'.format(spp_group))
        print("processing species group {}".format(spp_group))
        if not os.path.isfile(target_path):
            raster_list_sum(
                inputs_dir[spp_group], target_path, _TARGET_DATATYPE,
                nodata_remove=True)


def sum_terrestrial_mammals():
    """Sum up terrestrial mammals inputs for v2."""
    raster_list = [
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/felids_canids_sum_with_0.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data46/Bats_Distributions_Add.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data49/HaresPikas_Distributions_Add.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/land_mammals_with_caribou.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2/furbearers_sum_with_0.tif"]
    destination_folder = "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2"
    target_path = os.path.join(
        destination_folder, 'Terrestrial_Mammals_Habitats_add_v2.tif')
    raster_list_sum(
        raster_list, target_path, _TARGET_DATATYPE, nodata_remove=True)


def sum_terrestrial_inputs():
    """Sum up inputs to the terrestrial index."""
    destination_folder = "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2"
    raster_list = [
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data24/Amphibians_Habitat_v1.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data/AK_IBA.tif",
        os.path.join(destination_folder, 'Terrestrial_Mammals_Habitats_Final_v2.tif'),
        os.path.join(destination_folder, 'Sea_Birds_Habitat_v2.tif'),
        os.path.join(destination_folder, 'Land_Birds_Habitat_v2.tif')]
    target_path = os.path.join(
        destination_folder, 'AK_Terrestrial_Index_all_v2.tif')
    raster_list_sum(
        raster_list, target_path, _TARGET_DATATYPE, nodata_remove=True)


def sum_aquatic_inputs():
    """Sum up inputs to the aquatic index."""
    destination_folder = "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/species_groups_v2"
    raster_list = [
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/ShoreZone_intermediate/Marine_Habitats_final.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/efh_intermediate/EFH_not_in_SWAP_Final_v2.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/aquatic_index_revisions_8-9-21/Marine_Mammals_Final_v2.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data10/Lamprey_Habitats_Add.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data12/Pike_Char_Grayling_Habitat_Add.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data11/Whitefish_Blackfish_Inconnu_Reclass_Final.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data14/Rockfish_EFH_Final.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data15/Cod_Mackerel_Habitat_Final.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data17/Forage_Fish_Final_v1.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data16/Trout_Steelhead_Add.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data3/Salmon_Habitats_Final.tif",
        "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data8/Sablefish_Habitat_Final.tif",
        os.path.join(destination_folder, 'Crabs_Final_v2.tif'),
        os.path.join(destination_folder, 'marine_inverts_sum.tif')]
    target_path = os.path.join(
        destination_folder, 'AK_Aquatic_Index_all_v2.tif')
    raster_list_sum(
        raster_list, target_path, _TARGET_DATATYPE, nodata_remove=True)


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
        # 'Furbearers': [
        #     'AmericanBeaver_AnnualDistribution',
        #     'AmericanMarten_AnnualDistribution',
        #     'AmericanMink_AnnualDistribution',
        #     'CanadianLynx_AnnualDistribution',
        #     'Ermine_AnnualDistribution',
        #     'LeastWeasel_AnnualDistribution',
        #     'Muskrat_AnnualDistribution',
        #     'NorthAmericanPorcupine_AnnualDistribution',
        #     'RedFox_AnnualDistribution',
        #     'SnowshoeHare_AnnualDistribution',
        #     'Wolverine_AnnualDistribution',
        #     'Woodchuck_AnnualDistribution'],
        # 'Large_Land_Mammals': [
        #     'AmericanBlackBear_AnnualDistribution',
        #     'BrownBear_AnnualDistribution',
        #     'Coyote_AnnualDistribution',
        #     'DallsSheep_AnnualDistribution',
        #     'Moose_AnnualDistribution',
        #     'MountainGoat_AnnualDistribution',
        #     'Muskox_AnnualDistribution',
        #     'Wolf_AnnualDistribution'],
        }
    extract_folder = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Unzip"
    destination_folder = "E:/Datasets/AK_GAP_Analysis/Other_Subsistence_Spp/Processed"

    for sp_key in spp_group_dict:
        print("** working on species group {} **".format(sp_key))
        with tempfile.TemporaryDirectory(
                prefix=sp_key, dir=extract_folder) as temp_spp_dir:
            for sp in spp_group_dict[sp_key]:
                zip_path = os.path.join(_GAP_FOLDER, '{}.zip'.format(sp))
                dest_path = os.path.join(extract_folder, '{}.img'.format(sp))
                if os.path.isfile(dest_path):
                    os.remove(dest_path)
                with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_folder)

            # print("testing nodata values")
            # raster_list = [
            #     os.path.join(extract_folder, '{}.img'.format(sp)) for sp in
            #     spp_group_dict[sp_key]]
            # input_nodata = pygeoprocessing.get_raster_info(
            #     raster_list[0])['nodata'][0]
            # pixel_size = pygeoprocessing.get_raster_info(
            #     raster_list[0])['pixel_size']
            # for raster_path in raster_list:
            #     test_nodata = pygeoprocessing.get_raster_info(
            #         raster_path)['nodata'][0]
            #     if test_nodata != input_nodata:
            #         reclassify_nodata(raster_path, input_nodata)

            # print("aligning input rasters")
            # aligned_path_list = [
            #     os.path.join(temp_spp_dir, '{}.tif'.format(sp)) for sp in
            #     spp_group_dict[sp_key]]
            # pygeoprocessing.align_and_resize_raster_stack(
            #     raster_list, aligned_path_list, ['near'] * len(raster_list),
            #     pixel_size, 'union')

            # print("calculating sum")
            # intermediate_path = os.path.join(
            #     temp_spp_dir, '{}.tif'.format(sp_key))
            # raster_list_sum(
            #     aligned_path_list, input_nodata, intermediate_path,
            #     _TARGET_NODATA, _TARGET_DATATYPE, nodata_remove=True)

            # print("clipping to region")
            # target_path = os.path.join(
            #     destination_folder, '{}.tif'.format(sp_key))
            # pygeoprocessing.align_and_resize_raster_stack(
            #     [intermediate_path], [target_path], ['near'], pixel_size,
            #     'intersection', base_vector_path_list=[_BOUNDARY_PATH],
            #     raster_align_index=0,
            #     vector_mask_options={'mask_vector_path': _BOUNDARY_PATH})


def AK_terrestrial_index_sum():
    """Sum up inputs to the AK terrestrial index v2."""
    # TODO modify with species groups
    # TODO modify to use the new raster_list_sum()
    with tempfile.TemporaryDirectory() as temp_spp_dir:
        raster_list = [
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data/AK_IBA.tif",
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data24/Amphibians_Habitat_v1.tif",
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data25/Sea_Birds_Habitat_v1.tif",
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data26/Terrestrial_Mammals_Habitats_Final.tif",
            "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data22/Land_Birds_Habitat_v1.tif"]
        pixel_size = pygeoprocessing.get_raster_info(
            raster_list[0])['pixel_size']
        for raster_path in raster_list:
            test_nodata = pygeoprocessing.get_raster_info(
                raster_path)['nodata'][0]
            if test_nodata != _TARGET_NODATA:
                reclassify_nodata(raster_path, _TARGET_NODATA)

        print("aligning input rasters")
        aligned_path_list = [
            os.path.join(temp_spp_dir, os.path.basename(f)) for f in
            raster_list]
        pygeoprocessing.align_and_resize_raster_stack(
            raster_list, aligned_path_list, ['near'] * len(raster_list),
            pixel_size, 'union')

        print("calculating sum")
        intermediate_path = os.path.join(temp_spp_dir, 'terr_unclip.tif')
        raster_list_sum(
            aligned_path_list, _TARGET_NODATA, intermediate_path,
            _TARGET_NODATA, _TARGET_DATATYPE, nodata_remove=True)

        print("clipping to region")
        destination_path = "E:/Packages/AK_Wildlife_012721_4c7e75/commondata/raster_data26/AK_Terrestrial_Index_all_v2.tif"
        pygeoprocessing.align_and_resize_raster_stack(
            [intermediate_path], [destination_path], ['near'], pixel_size,
            'intersection', base_vector_path_list=[_BOUNDARY_PATH],
            raster_align_index=0,
            vector_mask_options={'mask_vector_path': _BOUNDARY_PATH})
        

if __name__ == "__main__":
    # gap_data_workflow()
    # add_elk_to_large_land_mammals()
    # AK_terrestrial_index_sum()
    # gap_data_workflow()
    # species_groups_sum()
    # sum_terrestrial_mammals()
    # sum_terrestrial_inputs()
    sum_aquatic_inputs()