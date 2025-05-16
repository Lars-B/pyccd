"""
WIP: Module to provide example code for the functions provided in the package.
"""
from pathlib import Path

from pyccd.read_nexus import read_nexus_trees


def read_transmission_nexus():
    """
    WIP read transmission nexus file
    :return:
    """
    # test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    test_tree_file = f"{Path(__file__).parent.absolute().parent}/tests/data/BREATH5taxa.trees"
    trees = read_nexus_trees(test_tree_file)
    # Testing how to output these trees, still WIP and unfinished
    # teststr = trees[0].write(
    #               features=["blockcount", "transm_ancest"], format_root_node=True, format=2)


def label_transmission_tree():
    """
    WIP label a given tree with its implied transmission tree
    :return:
    """
    tree_file = f"{Path(__file__).parent.absolute().parent}/tests/data/Filter-roetzer40.trees"
    trees = read_nexus_trees(tree_file)
    cur_tree_nwk = trees[1].write(features=["transm_ancest"], format_root_node=True, format=2)


if __name__ == '__main__':
    read_transmission_nexus()
    label_transmission_tree()
