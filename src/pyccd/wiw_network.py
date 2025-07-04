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
