from brokilon.ccd.domain.topology.ccd import get_clades
from brokilon.core.tree import TreeError, extract_edge_set


def robinson_foulds(tree_t1, tree_t2, attr_t1="name", attr_t2="name"):
    """
    Returns the Robinson-Foulds symmetric distance between current
    tree and a different tree instance.

    :param tree_t1: The first input tree for distance computation
    :param tree_t2: The second input tree for distance computation

    :param name attr_t1: Compare trees using a custom node
                          attribute as a node name.

    :param name attr_t2: Compare trees using a custom node
                          attribute as a node name in target tree.

    :returns: (rf, rf_max, common_attrs, names, edges_t1, edges_t2)
    """
    ref_t = tree_t1

    if len(ref_t.children) > 2 or len(tree_t2.children) > 2:
        raise TreeError(
            "Unrooted tree found! Not supported in this package...")

    attrs_t1 = {getattr(n, attr_t1) for n in ref_t.iter_leaves() if hasattr(n, attr_t1)}
    attrs_t2 = {getattr(n, attr_t2) for n in tree_t2.iter_leaves() if hasattr(n, attr_t2)}
    if attrs_t1 != attrs_t2:
        raise TreeError("Trees have different taxa sets, not supported...")

    edges1 = extract_edge_set(ref_t, attr_t1)
    edges2 = extract_edge_set(tree_t2, attr_t2)

    # the two root edges are never counted here, as they are always
    # present in both trees because of the common attr filters
    rf = len((edges1 ^ edges2))

    # Otherwise we need to count the actual number of valid
    # partitions in each tree -2 is to avoid counting the root
    # partition of the two trees (only needed in rooted trees)
    max_parts = (sum(1 for p in edges1 if len(p) > 1) +
                 sum(1 for p in edges2 if len(p) > 1) - 2)

    # min_comparison = (rf, max_parts, edges1, edges2)
    return rf


def robinson_foulds_ccd(tree_t1, tree_t2):
    # this is half as fast as the function above (which makes sense, its based on a weird function)
    clades_t1 = get_clades(tree_t1)
    clades_t2 = get_clades(tree_t2)
    return len(clades_t1.symmetric_difference(clades_t2))


if __name__ == '__main__':
    from pathlib import Path
    from brokilon.core import read_nexus_trees

    tree_file = (
        f"{Path(__file__).parent.absolute().parent.parent.parent}"
        f"/examples/data/400Taxa.trees"
    )
    trees = read_nexus_trees(tree_file)
    d1 = robinson_foulds(trees[0], trees[1], attr_t1="name", attr_t2="name")
    d2 = robinson_foulds_ccd(trees[0], trees[1])
    print(d1, d2)
