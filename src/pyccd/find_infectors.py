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

            # node_root_dist = node.get_distance(root_node)
            parent_root_dist = node.up.get_distance(root_node)
            # Using block end due to blocks, if no block this is the same as start
            event_length_from_parent = node.dist * node.blockend

            infection_start_dist_root = parent_root_dist + event_length_from_parent

            if node.blockcount > 0:
                block_name = f"block_{node.name}"
                data = [
                    block_name,
                    node.name,
                    infection_start_dist_root,
                    node.blockcount,
                ]
            else:
                data = [
                    node.transm_ancest,
                    node.name,
                    infection_start_dist_root,
                    None
                ]
            unknown_out_node = None
            if node.transm_ancest.startswith("Unknown"):
                unknown_out_node = node
            return node.transm_ancest, [data], unknown_out_node
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
            None,
        ]
        return parent.transm_ancest, [data], None
    # assert parent.blockstart == parent.blockend, "Assuming this is one infection atm!"
    # node_root_dist = node.get_distance(root_node)  # assumes node is a leaf!
    # blockend so that if its a block we have the correct start...
    length_to_add_from_cur_edge = parent.dist * parent.blockend
    infection_start_dist_root = parent.up.get_distance(root_node) + length_to_add_from_cur_edge
    if parent.blockcount > 0:
        block_name = f"block_{node.name}"
        data = [
            block_name,
            node.name,
            infection_start_dist_root,
            parent.blockcount,
        ]
    else:
        data = [
            parent.transm_ancest,
            node.name,
            infection_start_dist_root,
            None,
        ]
    unknown_out_node = None
    if parent.transm_ancest.startswith("Unknown"):
        unknown_out_node = parent
    return parent.transm_ancest, [data], unknown_out_node


def find_infector_unknown(cur_u_node, root_node):
    """
    Finds the infector of nodes that are unknown, not leafs

    :param cur_u_node: A node whose transm_ancest is unknown
    :param root_node: The root node of a tree (tree.get_tree_root())
    :return:
    """
    if cur_u_node.transm_ancest != cur_u_node.up.transm_ancest:
        if cur_u_node.up.is_root():
            parent_root_dist = 0.0
        else:
            parent_root_dist = cur_u_node.up.up.get_distance(root_node)

        event_start_length_from_parent = cur_u_node.up.dist * cur_u_node.up.blockend
        infection_start_dist_root = parent_root_dist + event_start_length_from_parent

        if cur_u_node.blockcount > 0:
            infection_start_dist_root = parent_root_dist + (cur_u_node.blockstart * cur_u_node.dist)
            block_name = f"block_{cur_u_node.name}"
            if cur_u_node.up.blockcount > 0:
                node_sibling = next(x for x in cur_u_node.up.children if x != cur_u_node)
                if node_sibling.blockcount > 0:
                    # infecting a block from a block with a block as sibling will result in a newly
                    # named unknown sample
                    data = [
                        # cur_u_node.up.transm_ancest,
                        f"Unknown-block_{cur_u_node.up.name}",
                        block_name,
                        infection_start_dist_root,
                        cur_u_node.blockcount,
                    ]
                else:
                    # infecting a block from a block with a sibling having some label we know that
                    # the sibling sample persists and infects this block
                    data = [
                        node_sibling.transm_ancest,
                        block_name,
                        infection_start_dist_root,
                        cur_u_node.blockcount,
                    ]
            else:
                data = [
                    cur_u_node.up.transm_ancest,
                    block_name,
                    infection_start_dist_root,
                    cur_u_node.blockcount,
                ]
        else:
            if cur_u_node.up.blockcount > 0:
                block_name = f"block_{cur_u_node.up.name}"
                data = [
                    block_name,
                    cur_u_node.transm_ancest,
                    infection_start_dist_root,
                    None,
                ]
            else:
                data = [
                    cur_u_node.up.transm_ancest,
                    cur_u_node.transm_ancest,
                    infection_start_dist_root,
                    None,
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
        data = [
            parent.transm_ancest,
            cur_u_node.transm_ancest,
            0.0,
            None,
        ]
        return [data]

    dist_start_infection_to_root = (parent.up.get_distance(root_node) +
                                    (parent.dist * parent.blockend))

    if cur_u_node.blockcount > 0:
        block_name = f"block_{cur_u_node.name}"
        if parent.blockcount > 0:
            node_sibling = next(x for x in cur_u_node.up.children if x != cur_u_node)
            if node_sibling.blockcount > 0:
                data = [
                    f"Unknown-block_{parent.name}",
                    block_name,
                    dist_start_infection_to_root,
                    cur_u_node.blockcount,
                ]
            else:
                data = [
                    node_sibling.transm_ancest,
                    block_name,
                    dist_start_infection_to_root,
                    cur_u_node.blockcount,
                ]
        else:
            data = [
                parent.transm_ancest,
                block_name,
                dist_start_infection_to_root,
                cur_u_node.blockcount
            ]
    else:
        if parent.blockcount > 0:
            block_name = f"block_{parent.name}"
            data = [
                block_name,
                cur_u_node.transm_ancest,
                dist_start_infection_to_root,
                None
            ]
        else:
            data = [
                parent.transm_ancest,
                cur_u_node.transm_ancest,
                dist_start_infection_to_root,
                None
            ]
    result = [data]
    if parent.transm_ancest.startswith("Unknown"):
        more_data = find_infector_unknown(parent, root_node)
        result.extend(more_data)
    return result
