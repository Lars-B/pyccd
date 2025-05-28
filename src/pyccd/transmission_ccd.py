"""
This module contains all the functions relevant for creating transmission CCDs.
"""
import random
import warnings
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum
from functools import total_ordering

import numpy as np

from .tree import Tree


class TypeCCD(Enum):
    """
    Enum representing different types of CCD.

    Attributes:
        BLOCKS: Represents the "Blocks" mode, whether there is a block event.
        ANCESTRY: Represents the "Ancestry" mode for processing transmission ancestry.
    """
    BLOCKS = "Blocks"
    ANCESTRY = "Ancestry"
    # able to add more in the future


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
class TransmissionBlockClade(BaseClade):
    """
    Clade with a flag indicating whether a transmission block (event) occurred on the edge.

    Attributes:
        has_block (bool): Whether this clade was preceded by a block of transmissions.
    """
    has_block: bool


@dataclass(frozen=True)
class TransmissionAncestryClade(BaseClade):
    """
    Clade with transmission ancestor, i.e. who infected this clade.

    Attributes:
        transm_ancest (str): Transmission ancestor
    """
    transm_ancest: str


def get_transmission_maps(trees: list[Tree] | tuple[Tree], type_str: str = "Ancestry") -> tuple:
    """
    Extracts all the relevant information from a list of Tree objects.
    The maps m1 and m2 are used as in the Larget approach for CCD1.
    With these we can construct a MAP tree, which can be annotated with
    branch lengths and blockcount summaries using the other two returns.

    :param type_str: Currently, either 'Blocks' or 'Ancestry' to determine the types of CCD to
                     construct.
    :param trees: list of Trees from which to extract the clade splits
                  from.
    :returns: Tuple of m1 (Clade counts), m2 (Clade split counts),
              blockcount_map (Blockcount counts), branch_lengths_map (Branch lengths)
    """
    try:
        # Converting type to an enum if possible
        ccd_type = TypeCCD(type_str)
    except ValueError as e:
        raise ValueError(f"Type '{type_str}' not recognized. "
                         f"Expected one of: {', '.join([item.value for item in TypeCCD])}") from e

    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(
        int)  # map for each (c1,c2) clade how often this specific relation got sampled

    blockcount_map = defaultdict(list)
    branch_lengths_map = defaultdict(list)

    # traversing all nodes of all trees using levelorder traversal
    for node in (node for t in trees for node in t.traverse("levelorder")):
        assert hasattr(node, "name"), "Node should have a name!"
        assert hasattr(node, "blockcount"), (
            "The nodes should have the blockcount "
            "attribute!")
        if len(node) > 1:
            assert len(node.children) == 2, "Non binary tree not supported!"

            parent_clade, child0_clade, child1_clade, blockcount_map, branch_lengths_map = (
                _add_internal_clade(node, ccd_type, blockcount_map,
                                    branch_lengths_map))

            m1[parent_clade] += 1
            # child0 and child1 are sorted within _add_internal_node function
            # child0 always contains the min leaf label among the two clades
            m2[(parent_clade, child0_clade, child1_clade)] += 1

        elif len(node) == 1:
            assert node.is_leaf(), "Should be a leaf node!"
            # leaf node for which we need to add the blockcount to the blockcount_map
            blockcount_map, branch_lengths_map = _add_leaf_clade(node, ccd_type,
                                                                 blockcount_map,
                                                                 branch_lengths_map)

    return m1, m2, blockcount_map, branch_lengths_map


def _add_internal_clade(node, ccd_type, blockcount_map: dict,
                        branch_lengths_map: dict) \
        -> tuple[BaseClade, BaseClade, BaseClade, dict, dict]:
    """
    Processes an internal node by constructing parent and child clades based on the CCD type,
    and updates blockcount and branch length maps accordingly.

    :param node: The internal tree node to process. Assumed to have two children, a blockcount,
                 branch length (dist), and optionally transmission ancestry.
    :param ccd_type: Type of CCD to use and construct clades for.
    :param blockcount_map: Dictionary mapping clades to a list of blockcounts.
    :param branch_lengths_map: Dictionary mapping clades to a list of branch lengths.
    :returns: A tuple containing the parent clade,
              the two child clades (ordered by minimum leaf label),
              and the updated blockcount and branch lengths maps.
    """
    c0_leafs = {int(leaf.name) for leaf in node.children[0]}
    c1_leafs = {int(leaf.name) for leaf in node.children[1]}
    parent_clade_set = frozenset(sorted(c0_leafs.union(c1_leafs)))

    match ccd_type:
        case TypeCCD.BLOCKS:
            parent_clade = TransmissionBlockClade(parent_clade_set,
                                                  node.blockcount != -1)
            child0_clade = TransmissionBlockClade(frozenset(c0_leafs),
                                                  node.children[
                                                      0].blockcount != -1)
            child1_clade = TransmissionBlockClade(frozenset(c1_leafs),
                                                  node.children[
                                                      1].blockcount != -1)
        case TypeCCD.ANCESTRY:
            parent_clade = TransmissionAncestryClade(parent_clade_set,
                                                     node.transm_ancest)
            child0_clade = TransmissionAncestryClade(frozenset(c0_leafs),
                                                     node.children[
                                                         0].transm_ancest)
            child1_clade = TransmissionAncestryClade(frozenset(c1_leafs),
                                                     node.children[
                                                         1].transm_ancest)
        case _:
            raise ValueError(f"Unknown type given: {ccd_type}")

    # Keeping track of blockcounts if not -1
    if node.blockcount != -1:
        blockcount_map[parent_clade].append(node.blockcount)

    # adding distance of parent clade to map of branch lengths
    branch_lengths_map[parent_clade].append(node.dist)

    # Depending on the lower value leaf we return child0 and child1 clades
    if min(c0_leafs) < min(c1_leafs):
        return parent_clade, child0_clade, child1_clade, blockcount_map, branch_lengths_map
    return parent_clade, child1_clade, child0_clade, blockcount_map, branch_lengths_map


def _add_leaf_clade(node, ccd_type: TypeCCD, blockcount_map: dict,
                    branch_lengths_map: dict) -> tuple[dict, dict]:
    """
    Processes a leaf node by creating an appropriate clade based on the CCD type,
    and updates blockcount and branch length maps accordingly.

    :param node: The tree node corresponding to a leaf.
    :param ccd_type: Type of CCD to use and construct clades for.
    :param blockcount_map: Dictionary mapping clades to a list of blockcounts.
    :param branch_lengths_map: Dictionary mapping clades to a list of branch lengths.
    :returns: Updated blockcount_map and branch_lengths_map.
    """
    match ccd_type:
        case TypeCCD.BLOCKS:
            leaf_clade = TransmissionBlockClade(frozenset({int(node.name)}),
                                                (node.blockcount != -1))
        case TypeCCD.ANCESTRY:
            leaf_clade = TransmissionAncestryClade(frozenset({int(node.name)}),
                                                   node.transm_ancest)
        case _:
            raise ValueError(f"Unknown type given: {ccd_type}")

    # Keeping track of blockcount for summarization, regardless of type
    if node.blockcount != -1:
        blockcount_map[leaf_clade].append(node.blockcount)
    # Adding the branch length to the branch_length_map
    branch_lengths_map[leaf_clade].append(node.dist)
    return blockcount_map, branch_lengths_map


def get_transmission_ccd_tree_bottom_up(m1: dict, m2: dict,
                                        blockcount_map: dict,
                                        branch_lengths_map: dict,
                                        seed: int = 42) -> str:
    """
    Constructs the transmission CCD MAP tree using a bottom-up approach.

    This function processes the given clades and splits, iterating over them
    in order of increasing clade size (smallest to largest). It calculates
    probabilities for each split, resolving ties randomly based on
    a 50% chance (controlled by a custom seed) and keeps only the highest
    probable splits for the MAP tree. The final resolved splits
    are used to construct the MAP tree, which is returned in Newick format.

    :param m1: A dictionary with clades as keys and their occurrences as values.
    :param m2: A dictionary with cladesplits as keys and their probabilities as values.
    :param blockcount_map: A mapping of clades to their respective block counts
    :param branch_lengths_map: A mapping of clades to their respective branch lengths
    :param seed: A seed for the random number generator to control tie-breaking. Default is 42.
    :returns: A string representing the tree in Newick format, annotated with median block-counts
              and mean branch lengths (might change in future versions).
    """
    # Seed used for tie breaking using random.random() <= 0.5
    random.seed(seed)

    seen_resolved_clades = {}

    # sorted list of all clades, small (cherries) to big
    for current_clade in sorted(list(m1.keys()), key=len):
        # the following are all triplets that represent how the current clade splits

        for current_split in (i for i in m2 if i[0] == current_clade):
            child1, child2 = current_split[1], current_split[2]

            assert current_clade == current_split[0], "This should be the same."
            assert (len(child1) == 1 or child1 in seen_resolved_clades), \
                "child1 should be in seen_resolved_clades when its length > 1"
            assert (len(child2) == 1 or child2 in seen_resolved_clades), \
                "child2 should be in seen_resolved_clades when its length > 1"

            c1_prob = 1 if len(child1) == 1 else seen_resolved_clades[child1][0]
            c2_prob = 1 if len(child2) == 1 else seen_resolved_clades[child2][0]

            # cur_prob = m2[current_split] / m1[current_split[0]]
            split_prob = c1_prob * c2_prob * (
                    m2[current_split] / m1[current_clade])

            if current_clade in seen_resolved_clades:
                if seen_resolved_clades[current_clade][0] < split_prob:
                    seen_resolved_clades[current_clade] = (split_prob,
                                                           current_split, False)
                elif seen_resolved_clades[current_clade][0] == split_prob:
                    # choose 50/50 if we want to update or keep the old better split.
                    if random.random() < 0.5:
                        # choose the new split
                        chosen_prob, chosen_split = split_prob, current_split
                    else:
                        # choose the old split
                        chosen_prob, chosen_split = seen_resolved_clades[current_clade][:2]
                    # Thrid entry True because tiebreaking is in effect for this split.
                    seen_resolved_clades[current_clade] = (chosen_prob, chosen_split, True)

            else:
                seen_resolved_clades[current_clade] = (split_prob,
                                                       current_split, False)

    # construct the root clade and build a tree dict
    root_clade = max(seen_resolved_clades.keys())
    output = _build_tree_dict_from_clade_splits(root_clade,
                                                seen_resolved_clades)

    return recursive_nwk_split_dict(root_clade, output, blockcount_map,
                                    branch_lengths_map)


def recursive_nwk_split_dict(clade, output, blockcount_map,
                             branch_lengths_map) -> str:
    """
    Recursively generates a Newick string for the given clade.
    Currently, it annotates the median blockcount if a block is present.
    If the given clade is a TransmissionAncestryClade it also annotates that.

    :param clade: The clade to generate the Newick string for.
    :param output: A dictionary containing the child clades for each parent.
                   As computed by the _build_tree_dirct_from_clade_splits function.
    :param blockcount_map: A dictionary mapping clades to their associated blockcount values.
    :param branch_lengths_map: A dictionary mapping clades to their branch lengths.
    :returns: A string representing the tree in Newick format,
              annotated with meadian blockcount and mean branch lengths.
    """
    if len(clade) == 1:
        # Base case for leaf node
        return (f"{next(iter(clade.clade))}"
                f"[&blockcount="
                f"{np.median(blockcount_map[clade]) if clade in blockcount_map else -1},"
                f"&transmission.ancestor="
                f"{clade.transm_ancest if isinstance(clade, TransmissionAncestryClade) else 'None'}"
                f"]"
                f":{np.mean(branch_lengths_map[clade])}")
    return (
        "("
        f"{recursive_nwk_split_dict(output[clade][0], output, blockcount_map, branch_lengths_map)},"
        f"{recursive_nwk_split_dict(output[clade][1], output, blockcount_map, branch_lengths_map)})"
        f"[&blockcount={np.median(blockcount_map[clade]) if clade in blockcount_map else -1},"
        f"&transmission.ancestor="
        f"{clade.transm_ancest if isinstance(clade, TransmissionAncestryClade) else 'None'}"
        f"]:"
        f"{np.mean(branch_lengths_map[clade])}"
    )


def _build_tree_dict_from_clade_splits(root_clade: BaseClade,
                                       seen_resolved_clades: dict) -> dict:
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
     :param seen_resolved_clades: A dictionary mapping clades (of any subclass of `BaseClade`)
                                  to a tuple:
                                  '(probability, (parent_clade, left_clade, right_clade))'.
                                  Only the split information is used for constructing the tree.
     :returns: A dictionary mapping each clade (of any subclass of `BaseClade`)
               to its child clades (left, right),
               which are also instances of a subclass of `BaseClade`.
     """
    stack = [root_clade]
    output = {}
    tiebreaking_occurred = False

    while stack:
        parent = stack.pop()
        _, (_, left, right), tiebreaking = seen_resolved_clades[parent]
        if tiebreaking:
            tiebreaking_occurred = True

        output[parent] = (left, right)

        if len(left) > 1:
            stack.append(left)
        if len(right) > 1:
            stack.append(right)

    # todo optional logging of which clades are affected by the tie break?
    if tiebreaking_occurred:
        warnings.warn("Tie breaking affected the constructed MAP tree!")
    return output
