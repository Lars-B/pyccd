from CCDpy.read_breath_nexus import read_transmission_nexus_history
from CCDpy.transmission_ccd import *
from pathlib import Path


# todo write proper tests for the functions...

def test_read_transmission_nexus():
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    trees = read_transmission_nexus_history(test_tree_file, False)
    # test = get_transmission_clades(trees[1], blockcountmap)
    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees)
    get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)
    assert False


def test_transmissionCCD_MAP_nexus():
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    output_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40-tCCD.tree"
    transmission_ccd_map_nexus(input_trees_file=test_tree_file,
                               output_tree_file=output_file,
                               overwrite=True,
                               burnin=.1)
    assert True


def test_lazy():
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/roetzer40.trees"
    output_file = f"{Path(__file__).parent.absolute()}/data/testing_tCCD.tree"
    transmission_ccd_map_nexus(input_trees_file=test_tree_file,
                               output_tree_file=output_file, overwrite=True, burnin=.1)
    assert True
