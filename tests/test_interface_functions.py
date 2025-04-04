"""
Testing interface functions
"""
import os
from pathlib import Path

from pyccd.interface_functions import transmission_ccd_map_nexus


def test_transmission_ccd_map_nexus():
    """
    Testing the transmission_ccd_map_nexus() interface function on a smaller (few taxa) tree set.
    """
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
    output_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40-tCCD.tree"
    transmission_ccd_map_nexus(input_trees_file=test_tree_file,
                               output_tree_file=output_file,
                               overwrite=True,
                               burn_in=.1)
    assert os.path.isfile(output_file)


def test_large_transmission_ccd_map_nexus():
    """
    Testing the transmission_ccd_map_nexus() interface function on a larger (more taxa) tree set.
    """
    test_tree_file = f"{Path(__file__).parent.absolute()}/data/roetzer40.trees"
    output_file = f"{Path(__file__).parent.absolute()}/data/testing_tCCD.tree"
    transmission_ccd_map_nexus(input_trees_file=test_tree_file,
                               output_tree_file=output_file, overwrite=True, burn_in=.1)
    assert os.path.isfile(output_file)
