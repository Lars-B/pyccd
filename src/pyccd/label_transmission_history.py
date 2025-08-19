"""
This module provides functions to label nodes in an Tree object
based on their blockcounts, assigning transmission ancestry identifiers.

It traverses the tree and assigns labels to nodes where blockcounts
indicate unknown transmission ancestry. Unlabeled nodes are tracked
for further processing.
"""
import collections
from typing import List, Tuple, Any

from .tree import Tree


def label_transmission_tree(tree):
    """
    Labels the transmission history onto a Tree object based on blockcounts.

    This function modifies the tree in place by assigning transmission ancestry
    labels to nodes where blockcounts indicate unknown transmission ancestry.
    It does not return any values.

    :param tree: A Tree object with blockcounts annotated.
    """

    unknown_count = 0

    unlabeled_nodes_list = []
    top_infected_nodes_list = collections.deque([])  # more efficient than list

    # Loops over all nodes first, could probably be done more efficiently in the future
    unlabeled_nodes_list, unknown_count = _label_all_nodes(
        tree, unlabeled_nodes_list, unknown_count
    )

    unlabeled_nodes_list, top_infected_nodes_list = (
        _label_leaf_and_reachable_nodes(tree, unlabeled_nodes_list, top_infected_nodes_list)
    )

    unlabeled_nodes_list, unknown_count = _label_top_infected_nodes(top_infected_nodes_list,
                                                                    unlabeled_nodes_list,
                                                                    unknown_count)

    if len(unlabeled_nodes_list) == 0:
        tree.get_tree_root().add_feature("num_unknowns", unknown_count)
        return

    # Sorting the list of unlabeled nodes by their level for _label_all_remaining_unknowns()
    unlabeled_nodes_list.sort(
        key=lambda obj: obj.get_distance(tree.get_tree_root(), topology_only=True))
    _label_all_remaining_unknowns(unlabeled_nodes_list, unknown_count)
    tree.get_tree_root().add_feature("num_unknowns", unknown_count)


def _label_leaf_and_reachable_nodes(tree, unlabeled_nodes_list, top_infected_nodes_list):
    """
        Labels all leaf nodes and recursively propagates the
        transmission ancestry to all reachable nodes.

        This function processes all the leaves in the provided tree (tree).
        For each leaf node:

        - If the blockcount is -1 (indicating no transmission event),
          it labels the leaf with its name and propagates this label upwards and to the
          children, marking all reachable nodes with the same ancestry.
        - Creates a queue of top infected nodes that arise during the propagation

        The function modifies the tree by adding the transmission ancestry feature to nodes and
        updates the `unlabeled_nodes_list` by removing nodes that have been labeled.

        :param tree: The tree to be processed (should contain leaf nodes with `blockcount`).
        :param unlabeled_nodes_list: List of nodes that have not yet been labeled
                                     with transmission ancestry.
        :param top_infected_nodes_list: List of nodes that are considered top-infected nodes,
                                        which are infected by transmission propagation.
        :returns: The updated `unlabeled_nodes_list` and `top_infected_nodes_list` after
                  propagation.
        """
    # First label all leafs and reachable nodes from leaves
    # todo too many nested blocks here, abstract and pull out into functions.
    for leaf in tree:
        if leaf.blockcount == -1:
            # No transmission event towards the current node,
            # recursively label all the nodes reachable from here
            leaf.transm_ancest = leaf.name
            if leaf in unlabeled_nodes_list:
                unlabeled_nodes_list.remove(leaf)
            working_node_list = [leaf.up]

            while working_node_list:
                working_node = working_node_list.pop()
                if working_node.blockcount == -1:  # and not working_node.is_root():
                    # There is no event towards current node
                    working_node.transm_ancest = leaf.name  # label the current node
                    if working_node in unlabeled_nodes_list:
                        unlabeled_nodes_list.remove(working_node)
                    # if the current node is not the root we can look further up
                    if working_node.up:
                        if not hasattr(working_node.up, "transm_ancest"):
                            # if the node upwards has not been labeled
                            # we can check if we need to add it
                            if working_node.up not in working_node_list:
                                working_node_list.append(working_node.up)
                # If working_node has children we need to sort those out too
                if len(working_node.children) > 0:
                    assert len(working_node.children) == 2, "Non binary tree found, not supported!"
                    for c in working_node.children:
                        # adding the children that can be reached from current node
                        if not hasattr(c, "transm_ancest"):
                            if c.blockcount == -1:
                                # adding the children that are reachable
                                # but have not been labeled yet
                                if c not in working_node_list:
                                    working_node_list.append(c)
                            elif c.blockcount == 0:
                                # child that is infected by current label
                                c.transm_ancest = leaf.name
                                top_infected_nodes_list.append(c)
                                if c in unlabeled_nodes_list:
                                    # Removing the child c from unlabeled list
                                    unlabeled_nodes_list.remove(c)
        elif (leaf.blockcount == 0 and
              leaf not in top_infected_nodes_list and
              leaf not in unlabeled_nodes_list):
            unlabeled_nodes_list.append(leaf)
    return unlabeled_nodes_list, top_infected_nodes_list


def _label_all_nodes(tree: Tree, unlabeled_nodes_list: List = None,
                     unknown_count: int = 0) \
        -> \
                Tuple[List, int]:
    """
    Labels all nodes by iterating over the entire tree.

    :param tree: The tree to label.
    :param unlabeled_nodes_list: List of unlabeled nodes. If None, initializes to an empty list.
    :param unknown_count: Number of labeled unknown nodes (unknown transmission ancestors).
    :returns: A tuple containing the updated list of unlabeled nodes and the updated unknown count.
    """

    if unlabeled_nodes_list is None:
        unlabeled_nodes_list = []

    for node in tree.traverse("levelorder"):
        # Skip nodes that already have a transmission ancestor
        if hasattr(node, "transm_ancest"):
            continue

        if node.is_root():
            # Enforce root blockcount to be -1 (ensure this is correct)
            node.blockcount = -1

        if node.blockcount > 0:
            # Node has a block, assign an unknown transmission ancestor label
            node.transm_ancest = f"Unknown-{unknown_count}"
            unknown_count += 1

        else:
            unlabeled_nodes_list.append(node)

    return unlabeled_nodes_list, unknown_count


def _label_top_infected_nodes(top_infected_nodes_list, unlabeled_nodes_list, unknown_count) -> \
tuple[Any, Any]:
    """
    Labels transmission ancestry for top-infected nodes and their children.

    This function iterates over the list of top-infected nodes and propagates
    the transmission ancestry down to their children. If a child node is
    unlabeled and has a blockcount of 0 or -1, it inherits the transmission
    ancestry from its parent (the current node). If the child has a blockcount
    of -1, it is added to the list of top-infected nodes to process further.

    The function ensures that:
    - Each node in the `top_infected_nodes_list` is properly labeled.
    - Non-binary trees are not processed (only binary trees are supported).
    - Child nodes are removed from the `unlabeled_nodes_list` once labeled.

    :param top_infected_nodes_list: List of nodes with transmission ancestry
                                    that need to propagate labels to their children.
    :param unlabeled_nodes_list: List of nodes that need to be labeled with
                                 transmission ancestry.
    :returns: Updated `unlabeled_nodes_list` after labeling top-infected nodes' children.
    :raises AssertionError: If an unlabeled node is found or a non-binary tree is encountered.
    """
    while top_infected_nodes_list:
        cur_node = top_infected_nodes_list.pop()
        assert hasattr(cur_node, "transm_ancest"), "This node has to be labeled at this point.!"
        assert cur_node not in unlabeled_nodes_list, \
            "This node should not be in the unlabeled list!"
        if len(cur_node.children) > 0:
            added_unknown = False
            assert len(cur_node.children) == 2, "Non binary tree found, not supported!"
            child1, child2 = cur_node.children
            # pull down label of parent if children are reachable and unlabeled
            if not hasattr(child1, "transm_ancest"):
                if child1.blockcount == -1:
                    child1.transm_ancest = cur_node.transm_ancest
                    top_infected_nodes_list.append(child1)
            if not hasattr(child2, "transm_ancest"):
                if child2.blockcount == -1:
                    child2.transm_ancest = cur_node.transm_ancest
                    top_infected_nodes_list.append(child2)

            if not hasattr(child1, "transm_ancest"):
                if not hasattr(child2, "transm_ancest"):
                    # both children are unlabeld and have a blockcount > -1
                    # this is an unknown transmission ancestor
                    assert child1.blockcount > -1, ("Has to have at least "
                                                    "one transmission even on branch")
                    assert child2.blockcount > -1, ("Has to have at least "
                                                    "one transmission even on branch")
                    child1.transm_ancest = f"Unknown-{unknown_count}"
                    child2.transm_ancest = f"Unknown-{unknown_count}"
                    unknown_count += 1
                # only child1 is unlabeled and the infector of child2 is the infector of child1
                child1.transm_ancest = child2.transm_ancest
            if not hasattr(child2, "transm_ancest"):
                # only child 2 is still unlabled
                child2.transm_ancest = child1.transm_ancest

            if child1 in unlabeled_nodes_list:
                unlabeled_nodes_list.remove(child1)
            if child2 in unlabeled_nodes_list:
                unlabeled_nodes_list.remove(child2)
    return unlabeled_nodes_list, unknown_count


def _label_all_remaining_unknowns(unlabeled_nodes_list: List, unknown_count: int) -> None:
    """
    Labels all nodes in the provided list with the unknown transmission history
    IMPORTANT: This function assumes that the list of unlabeled nodes provided
    is sorted according to their level in the tree.

    The function modifies the nodes in place and does not return
    any values. The unknown_count is updated as each node is labeled.

    The method works as follows:

    - If a node is the root: it is labeled with the current unknown_count.
    - If a node's parent (up) has a blockcount of -1:
      the node inherits the parent’s transmission ancestry.
    - If neither of the above conditions are true:
      the function checks the node’s sibling for a label.
    - If the sibling has a transmission ancestry:
      the node is labeled with the same identifier.
    - If no valid label is found:
      the node is labeled with a new "Unknown-{unknown_count}" label.

    :param unlabeled_nodes_list: List of nodes to be labeled
                                 with transmission ancestry identifiers.
    :param unknown_count: The current count of labeled "Unknown" nodes
    """
    for node in unlabeled_nodes_list:
        if hasattr(node, "transm_ancest"):
            # This case should not happen, if it does it needs to be figured out...
            raise NotImplementedError("Needs correct implementation... or bug fix...")
        assert node.blockcount < 1, "These nodes should be labeled at this point in the code."
        # all the paths from leaves have been labeled
        # at this point there is no path to any leaf from here
        # since we have ordered the nodes in the list
        # according to their 'level' in the tree the ancestry pulls down
        if node.is_root():
            node.transm_ancest = f"Unknown-{unknown_count}"
            unknown_count += 1
        elif node.up.blockcount == -1:
            node.transm_ancest = node.up.transm_ancest
        else:
            node_sibling = next(x for x in node.up.children if x != node)
            if node_sibling.blockcount == -1:
                if hasattr(node_sibling, "transm_ancest"):
                    node.transm_ancest = node_sibling.transm_ancest
                else:
                    # Internal path that consists of unknown transmission
                    # ancestry but haven't labeled it yet
                    node.transm_ancest = f"Unknown-{unknown_count}"
                    unknown_count += 1
            if node_sibling.blockcount == 0 and hasattr(node_sibling, "transm_ancest"):
                # Path of unknown transmission, if already labeled reuse that unknown label!
                node.transm_ancest = node_sibling.transm_ancest
            else:
                node.transm_ancest = f"Unknown-{unknown_count}"
                unknown_count += 1
