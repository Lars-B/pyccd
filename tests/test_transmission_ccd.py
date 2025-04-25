"""
Tests for the transmission_ccd.py functions
"""
from pathlib import Path

import pytest

import pyccd.read_nexus
from pyccd.transmission_ccd import (TypeCCD, get_transmission_maps,
                                    get_transmission_ccd_tree_bottom_up)


def test_valid_enum():
    """
    Testing valid enum creation for TypeCCD
    """
    assert TypeCCD("Blocks") == TypeCCD.BLOCKS
    assert TypeCCD("Ancestry") == TypeCCD.ANCESTRY


def test_invalid_enum():
    """
    Testing invalid enum creation for TypeCCD
    """
    with pytest.raises(ValueError):
        TypeCCD("blocks")


def test_get_transmission_maps_wrong_enum_conversion():
    """
    Testing invalid type_str input for get_transmission_maps()
    """
    with pytest.raises(ValueError) as e:
        get_transmission_maps(trees=[], type_str="Bla")
    assert str(e.value).startswith("Type 'Bla' not recognized.")


def test_read_transmission_nexus():
    """
        WIP This will be moved to example folder shortly
    """
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    trees = pyccd.read_nexus.read_nexus_trees(test_tree_file, breath_trees=True,
                                              label_transm_history=False)
    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees)
    get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)
    assert True, "Just testing if nothing breaks, apparently something did..."


def test_transmisison_ccd_hiostory():
    """
        WIP This will be moved to example folder shortly
    """
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    trees = pyccd.read_nexus.read_nexus_trees(test_tree_file, breath_trees=True,
                                              label_transm_history=True)
    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees, type_str="Ancestry")
    tree = get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)
    assert tree is not None, "Failed this example..."
