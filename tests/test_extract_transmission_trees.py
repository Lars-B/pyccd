from CCDpy.extract_transmission_trees import *
from pathlib import Path


def test_read_transmission_nexus():
    # test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/BREATH5taxa.trees"
    trees = read_transmission_nexus(test_tree_file)
    # Testing how to output these trees, still WIP and unfinished
    # teststr = trees[0].write(features=["blockcount", "transm_ancest"], format_root_node=True, format=2)
    assert False


def test_label_transmission_tree():
    tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    trees = read_transmission_nexus(tree_file)
    assert False
