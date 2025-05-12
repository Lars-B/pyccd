"""
WIP.
Todo make this a CLI tool ?!?
This module contains easy to use interface functions, minimal input to get a desired output.
"""
import os

from pyccd.read_nexus import read_nexus_trees
from pyccd.transmission_ccd import get_transmission_ccd_tree_bottom_up, get_transmission_maps


def transmission_ccd_map_nexus(input_trees_file: str, output_tree_file: str,
                               overwrite: bool = False, burn_in: float = 0):
    """
    Todo: weird because it does everything in one go, maybe rename and use in conjuction with all
     possible ccds as an interface just to get the MAP tree?

    Takes input trees file and writes a tCCD1-MAP tree to the output file in Nexus format.

    :param input_trees_file: Path to the input trees file
    :type input_trees_file: str
    :param output_tree_file: Path to the output trees file
    :type output_tree_file: str
    :param overwrite: If true, overwrite existing tCCD-MAP trees
    :type overwrite: bool
    :param burn_in: Value between 0 and 1, defining how many trees to discard as burn-in
    :type burn_in: float
    :return: None
    """
    if os.path.exists(output_tree_file):
        if not overwrite:
            raise FileExistsError("File exists, set 'overwrite=True' to enable overwriting.")
        os.remove(output_tree_file)

    if not os.path.exists(input_trees_file):
        raise FileNotFoundError(f"Input trees {input_trees_file} not found.")

    if not 0 <= burn_in < 1:
        raise ValueError("Burnin should be a number between 0 and 1, representing the proportion "
                         "of trees to delete at the beginning of the file.")

    trees = read_nexus_trees(input_trees_file, breath_trees=True)
    trees = trees[int(burn_in * len(trees)):]  # deleting burn_in from trees
    if len(trees) < 1:
        raise ValueError("Treeset is empty after burn_in removal... reduce burning or check file.")

    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees)
    newick_map = get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)

    with (open(input_trees_file, 'r', encoding="UTF-8") as infile,
          open(output_tree_file, 'w+', encoding="UTF-8") as outfile):
        for line in infile:
            if line.strip().startswith("tree "):
                break  # Stop reading when the\ tree section starts
            outfile.write(line)

        outfile.write(f"tree tCCD_MAP = {newick_map};\nEnd;\n")
