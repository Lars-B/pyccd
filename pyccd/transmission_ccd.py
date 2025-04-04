"""
This module contains all the functions relevant for creating transmission CCDs.
"""
import os.path
import random
import warnings
from collections import defaultdict
from dataclasses import dataclass
from functools import total_ordering

import numpy as np

from pyccd.read_nexus import read_nexus_trees


@total_ordering
@dataclass(frozen=True)
class BaseClade:
    """
    Base class representing a clade — a set of taxa/leaves.

    Attributes:
        clade (frozenset): A frozen set of taxa or node labels in this clade.
    """
    clade: frozenset

    def __len__(self) -> int:
        """Return the number of taxa/nodes in the clade."""
        return len(self.clade)

    def __lt__(self, other) -> bool:
        """Compare clade sizes for sorting or ordering."""
        if isinstance(other, BaseClade):
            return len(self.clade) < len(other.clade)
        return NotImplemented

    def __eq__(self, other) -> bool:
        """Clades are equal if their taxon sets are equal."""
        return isinstance(other, BaseClade) and self.clade == other.clade


@dataclass(frozen=True)
class TransmissionClade(BaseClade):
    """
    Clade with a flag indicating whether a transmission block (event) occurred on the edge.

    Attributes:
        has_block (bool): Whether this clade was preceded by a block of transmissions.
    """
    has_block: bool


# @dataclass(frozen=True)
# class SupportClade(BaseClade):
#     """
#     Clade with a support value (e.g., posterior or bootstrap probability).
#
#     Attributes:
#         support (float): Confidence value assigned to this clade.
#     """
#     support: float


def get_transmission_maps(trees: list) -> tuple:
    """
    Extracts all the relevant information from a list of ete3.Tree objects.
    The maps m1 and m2 are used as in the Larget approach for CCD1.
    With these we can construct a MAP tree, which can be annotated with
    branch lengths and blockcount summaries using the other two returns.

    :param trees: list of ete3.Trees from which to extract the clade splits from.
    :type trees: list

    :returns: Tuple of m1 (Clade counts), m2 (Clade split counts),
              blockcount_map (Blockcount counts), branch_lengths_map (Branch lengths)
    :rtype: tuple
    """
    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled

    blockcount_map = defaultdict(list)
    branch_lengths_map = defaultdict(list)

    # traversing all nodes of all trees using levelorder traversal
    for node in (node for t in trees for node in t.traverse("levelorder")):
        if len(node) > 1:
            assert len(node.children) == 2, "Non binary tree not supported!"
            assert hasattr(node, "name"), "Node should have a name!"
            assert hasattr(node, "blockcount"), ("The nodes should have the blockcount "
                                                 "attribute!")

            c0_leafs = {int(leaf.name) for leaf in node.children[0]}
            c1_leafs = {int(leaf.name) for leaf in node.children[1]}
            parent_clade_set = frozenset(sorted(c0_leafs.union(c1_leafs)))

            if node.blockcount == -1:
                parent_clade = TransmissionClade(parent_clade_set, False)
            else:
                parent_clade = TransmissionClade(parent_clade_set, True)
                blockcount_map[parent_clade].append(node.blockcount)

            # adding distance of parent clade to map of branch lengths
            branch_lengths_map[parent_clade].append(node.dist)

            m1[parent_clade] += 1

            child0_clade = TransmissionClade(frozenset(c0_leafs), node.children[0].blockcount != -1)
            child1_clade = TransmissionClade(frozenset(c1_leafs), node.children[1].blockcount != -1)

            # in m2 the split clade with the lower int for taxa is entered first
            if min(c0_leafs) < min(c1_leafs):
                m2[(parent_clade, child0_clade, child1_clade)] += 1
            else:
                m2[(parent_clade, child1_clade, child0_clade)] += 1
        elif len(node) == 1:
            assert node.is_leaf(), "Should be a leaf node!"
            # leaf node for which we need to add the blockcount to the blockcount_map

            leaf_clade = TransmissionClade(frozenset({int(node.name)}), (node.blockcount != -1))
            # adding branch length of the leaf node to the dict
            branch_lengths_map[leaf_clade].append(node.dist)

            # if the clade has a block we need to keep track of the blockcount for summaries
            if node.blockcount != -1:
                blockcount_map[leaf_clade].append(node.blockcount)

    return m1, m2, blockcount_map, branch_lengths_map


def get_transmission_ccd_tree_bottom_up(m1: dict, m2: dict, blockcount_map: dict,
                                        branch_lengths_map: dict, seed: int = 42) -> str:
    """
    Constructs the transmission CCD MAP tree using a bottom-up approach.

    This function processes the given clades and splits, iterating over them
    in order of increasing clade size (smallest to largest). It calculates
    probabilities for each split, resolving ties randomly based on
    a 50% chance (controlled by a custom seed) and keeps only the highest
    probable splits for the MAP tree. The final resolved splits
    are used to construct the MAP tree, which is returned in Newick format.

    :param m1: A dictionary with clades as keys and their occurrences as values.
    :type m1: dict
    :param m2: A dictionary with cladesplits as keys and their probabilities as values.
    :type m2: dict
    :param blockcount_map: A mapping of clades to their respective block counts
    :type blockcount_map: dict
    :param branch_lengths_map: A mapping of clades to their respective branch lengths
    :type branch_lengths_map: dict
    :param seed: A seed for the random number generator to control tie-breaking. Default is 42.
    :type seed: int

    :returns: A string representing the tree in Newick format, annotated with median blockcounts
              and mean branch lengths (might change in future versions).
    :rtype: str
    """
    # Seed used for tie breaking using random.random() <= 0.5
    random.seed(seed)

    seen_resolved_clades = {}

    # sorted list of all clades, small (cherries) to big
    for current_clade in sorted(list(m1.keys()), key=len):
        # the following are all triplets that represent how the current clade splits

        for current_split in (i for i in m2 if i[0] == current_clade):
            child1, child2 = current_split[1], current_split[2]

            assert (len(child1) == 1 or child1 in seen_resolved_clades), \
                "child1 should be in seen_resolved_clades when its length > 1"
            assert (len(child2) == 1 or child2 in seen_resolved_clades), \
                "child2 should be in seen_resolved_clades when its length > 1"

            c1_prob = 1 if len(child1) == 1 else seen_resolved_clades[child1][0]
            c2_prob = 1 if len(child2) == 1 else seen_resolved_clades[child2][0]

            # cur_prob = m2[current_split] / m1[current_split[0]]
            split_prob = c1_prob * c2_prob * (m2[current_split] / m1[current_split[0]])

            if current_split[0] in seen_resolved_clades:
                if seen_resolved_clades[current_split[0]][0] < split_prob:
                    seen_resolved_clades[current_split[0]] = (split_prob, current_split)
                elif seen_resolved_clades[current_split[0]][0] == split_prob:
                    warnings.warn("Tie breaking in effect.")  # currently warns about tie breaking
                    if random.random() < 0.5:
                        # choose 50/50 if we want to update or keep the old better split.
                        seen_resolved_clades[current_split[0]] = (split_prob, current_split)
            else:
                seen_resolved_clades[current_split[0]] = (split_prob, current_split)

    # construct the root clade and build a tree dict
    root_clade = max(seen_resolved_clades.keys())
    output = _build_tree_dict_from_clade_splits(root_clade, seen_resolved_clades)

    return recursive_nwk_split_dict(root_clade, output, blockcount_map, branch_lengths_map)


def recursive_nwk_split_dict(clade, output, blockcount_map, branch_lengths_map):
    """
    Recursively generates a Newick string for the given clade.
    Currently, it annotates the median blockcount if a block is present.

    :param clade: The clade to generate the Newick string for.
    :type clade: TransmissionClade
    :param output: A dictionary containing the child clades for each parent.
                   As computed by the _build_tree_dirct_from_clade_splits function.
    :type output: dict
    :param blockcount_map: A dictionary mapping clades to their associated blockcount values.
    :type blockcount_map: dict
    :param branch_lengths_map: A dictionary mapping clades to their branch lengths.
    :type branch_lengths_map: dict
    :return: A string representing the tree in Newick format,
             annoated with meadian blockcount and mean branch lengths.
    :rtype: str
    """
    if len(clade) == 1:
        # Base case for leaf node
        return (f"{next(iter(clade.clade))}"
                f"[&blockcount={np.median(blockcount_map[clade]) if clade.has_block else -1}]"
                f":{np.mean(branch_lengths_map[clade])}")
    # recursive case for internal node
    return (f"({recursive_nwk_split_dict(output[clade][0], output,
                                         blockcount_map, branch_lengths_map)},"
            f"{recursive_nwk_split_dict(output[clade][1], output,
                                        blockcount_map, branch_lengths_map)})"
            f"[&blockcount={np.median(blockcount_map[clade]) if clade.has_block else -1}]"
            f":{np.mean(branch_lengths_map[clade])}")


def _build_tree_dict_from_clade_splits(root_clade: BaseClade, seen_resolved_clades: dict) -> dict:
    """
     Constructs a tree dictionary from a set of resolved clade splits.

     Given a root clade and a dictionary of previously resolved clade splits,
     this function recursively builds a dictionary representing the binary
     tree structure. Each entry maps a parent clade (of any subclass of `BaseClade`)
     to a tuple of its left and right child clades.

     This function is designed to work with any class that inherits from the base
     `BaseClade` class, making it flexible for future extensions to other clade types.

     :param root_clade: The root clade to start building the tree from.
                        This clade and its children are instances of a subclass of `Clade`.
     :type root_clade: BaseClade
     :param seen_resolved_clades: A dictionary mapping clades (of any subclass of `BaseClade`)
                                  to a tuple:
                                  '(probability, (parent_clade, left_clade, right_clade))'.
                                  Only the split information is used for constructing the tree.
     :type seen_resolved_clades: dict
     :return: A dictionary mapping each clade (of any subclass of `BaseClade`)
              to its child clades (left, right),
              which are also instances of a subclass of `BaseClade`.
     :rtype: dict
     """
    stack = [root_clade]
    output = {}

    while stack:
        parent = stack.pop()
        _, (_, left, right) = seen_resolved_clades[parent]
        # cur_split = seen_resolved_clades[cur_parent][1]
        output[parent] = (left, right)
        # output[cur_parent] = (cur_split[1], cur_split[2])
        if len(left) > 1:
            stack.append(left)
        if len(right) > 1:
            stack.append(right)
    return output


def transmission_ccd_map_nexus(input_trees_file: str, output_tree_file: str,
                               overwrite: bool = False, burnIn: float = 0):
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
    :param burnIn: Value between 0 and 1, defining how many trees to discard as burn-in
    :type burnIn: float
    :return: None
    """
    if os.path.exists(output_tree_file):
        if not overwrite:
            raise FileExistsError("File exists, set 'overwrite=True' to enable overwriting.")
        os.remove(output_tree_file)

    if not os.path.exists(input_trees_file):
        raise FileNotFoundError(f"Input trees {input_trees_file} not found.")

    if not 0 <= burnIn < 1:
        raise ValueError("Burnin should be a number between 0 and 1, representing the proportion "
                         "of trees to delete at the beginning of the file.")

    trees = read_nexus_trees(input_trees_file, breath_trees=True)
    trees = trees[int(burnIn * len(trees)):]  # deleting burnIn from trees
    if len(trees) < 1:
        raise ValueError("Treeset is empty after burnIn removal... reduce burning or check file.")

    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees)
    newick_map = get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)

    with (open(input_trees_file, 'r', encoding="UTF-8") as infile,
          open(output_tree_file, 'w+', encoding="UTF-8") as outfile):
        for line in infile:
            if line.strip().startswith("tree "):
                break  # Stop reading when the\ tree section starts
            outfile.write(line)

        outfile.write(f"tree tCCD_MAP = {newick_map};\nEnd;\n")
