"""
Module contains functions regarding the computation and plotting of WIW networks for transmission
trees and transmission CCDs.
"""


def find_infector(node, indirect: bool = False):
    """
    Takes a node in a transmission annotated phylogenetic tree and returns its infector.

    :param indirect: Infer indirect transmissions, i.e. find last known infector
    :param node: Node of which to find the infector.
    :return: Infector of node
    """
    if node.name != node.transm_ancest:
        if not (indirect and node.transm_ancest.startswith("Unknown")):
            return node.transm_ancest
    parent = node.up
    while (parent.transm_ancest == node.name or
           (indirect and parent.transm_ancest.startswith("Unknown"))):
        if parent.is_root():
            break
        parent = parent.up
    if indirect and (not parent.is_root()) and parent.transm_ancest.startswith("Unknown"):
        raise ValueError("This should not happen?")
    return parent.transm_ancest


def find_infector_with_data(node, root_node, indirect: bool = False):
    """
    Finds the infector with extra data: distance to root of start and end of infection

    :param node: The node to find the infector of, in this case assumes to be a leaf
    :param root_node: The root node object of a tree (tree.get_tree_root())
    :param indirect: Currently not supported?
    :return:
    """
    if node.name != node.transm_ancest:
        if not (indirect and node.transm_ancest.startswith("Unknown")):
            # regular infection from one taxon to the next or Unknown/Block

            node_root_dist = node.get_distance(root_node)
            parent_root_dist = node.up.get_distance(root_node)
            assert node_root_dist > parent_root_dist, "Impossible!"
            # Using block end due to blocks, if no block this is the same as start
            event_length_from_parent = node.dist * node.blockend

            infection_start_dist_root = parent_root_dist + event_length_from_parent

            data = [
                node.transm_ancest,
                node.name,
                infection_start_dist_root,
                node_root_dist,
                node.blockcount
            ]
            unknown_out_node = None
            if node.transm_ancest.startswith("Unknown"):
                unknown_out_node = node
            return node.transm_ancest, data, unknown_out_node
    parent = node.up
    while (parent.transm_ancest == node.name or
           (indirect and parent.transm_ancest.startswith("Unknown"))):
        if parent.is_root():
            break
        parent = parent.up
    if indirect and (not parent.is_root()) and parent.transm_ancest.startswith("Unknown"):
        raise ValueError("This should not happen?")
    if parent.is_root():
        # parent is the root is the leaf, og infection
        assert node.name == parent.transm_ancest, "The parent is the root and node.name is OG?"
        assert node.blockcount == -1, "Something wrong?"
        data = [
            parent.transm_ancest,
            node.name,
            0.0,
            node.get_distance(root_node),
            node.blockcount
        ]
        return parent.transm_ancest, data, None
    # assert parent.blockstart == parent.blockend, "Assuming this is one infection atm!"
    node_root_dist = node.get_distance(root_node)  # assumes node is a leaf!
    # blockend so that if its a block we have the correct start...
    length_to_add_from_cur_edge = parent.dist * parent.blockend
    infection_start_dist_root = parent.up.get_distance(root_node) + length_to_add_from_cur_edge
    data = [
        parent.transm_ancest,
        node.name,
        infection_start_dist_root,
        node_root_dist,
        node.blockcount
    ]
    unknown_out_node = None
    if parent.transm_ancest.startswith("Unknown"):
        unknown_out_node = parent
    return parent.transm_ancest, data, unknown_out_node


def find_infector_unknown(cur_u_node, root_node):
    """
    Finds the infector of nodes that are unknown, not leafs

    :param cur_u_node: A node whoese transm_ancest is unknown
    :param root_node: The root node of a tree (tree.get_tree_root())
    :return:
    """
    if cur_u_node.transm_ancest != cur_u_node.up.transm_ancest:
        parent_root_dist = cur_u_node.up.get_distance(root_node)
        event_start_length_from_parent = cur_u_node.dist * cur_u_node.blockstart
        event_end_length_from_parent = cur_u_node.dist * cur_u_node.blockend

        infection_start_dist_root = parent_root_dist + event_start_length_from_parent
        infection_end_dist_root = parent_root_dist + event_end_length_from_parent

        data = [
            cur_u_node.up.transm_ancest,
            cur_u_node.transm_ancest,
            infection_start_dist_root,
            infection_end_dist_root,
            cur_u_node.blockcount
        ]
        result = [data]
        if cur_u_node.up.transm_ancest.startswith("Unknown"):
            if not cur_u_node.up.is_root():
                more_data = find_infector_unknown(cur_u_node.up, root_node)
                result.extend(more_data)
        return result
    parent = cur_u_node.up
    while parent.transm_ancest == cur_u_node.transm_ancest:
        if parent.is_root():
            break
        parent = parent.up
    if parent.is_root():
        assert parent.transm_ancest == cur_u_node.transm_ancest, "This should not happen!"
        parent_root_dist = cur_u_node.up.get_distance(root_node)
        event_end_length_from_parent = cur_u_node.dist * cur_u_node.blockstart
        data = [
            parent.transm_ancest,
            cur_u_node.transm_ancest,
            0.0,
            parent_root_dist + event_end_length_from_parent,
            cur_u_node.blockcount
        ]
        return [data]
    dist_end_infection_to_root = (cur_u_node.up.get_distance(root_node) +
                                  (cur_u_node.dist * cur_u_node.blockstart))

    dist_start_infection_to_root = (parent.up.get_distance(root_node) +
                                    (parent.dist * parent.blockend))

    data = [
        parent.transm_ancest,
        cur_u_node.transm_ancest,
        dist_start_infection_to_root,
        dist_end_infection_to_root,
        cur_u_node.blockcount
    ]
    result = [data]
    if parent.transm_ancest.startswith("Unknown"):
        more_data = find_infector_unknown(parent, root_node)
        result.extend(more_data)
    return result
