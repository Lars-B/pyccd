from CCDpy.transmissionCCD import *
from pathlib import Path


def test_read_transmission_nexus():
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    trees = read_transmission_nexus(test_tree_file)
    # test = get_transmission_clades(trees[1], blockcountmap)
    m1, m2, blockcount_map = get_transmission_maps(trees)
    get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map)
    assert False


def test_transmissionCCD_MAP_nexus():
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    output_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40-tCCD.tree"
    transmissionCCD_MAP_nexus(input_trees_file=test_tree_file,
                              output_tree_file=output_file, overwrite=True)
    assert True
