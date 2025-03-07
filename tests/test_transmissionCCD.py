from CCDpy.transmissionCCD import *
from pathlib import Path


def test_read_transmission_nexus():
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/BREATH5taxa.trees"
    trees = read_transmission_nexus(test_tree_file)
    blockcountmap = {}
    test = get_transmission_clades(trees[1], blockcountmap)
    assert False


def test_transmission_clades():
    test = transmission_clade(frozenset({1, 2}), False)
    test1 = transmission_clade(frozenset({1, 2}), True)
    testingdict = defaultdict()
    testingdict[test] = 1
    testingdict[test1] = 5
    # todo this works, implement the same for the other function
    #  need to add dictionary to keep track of the values for block count probably?
    assert False
