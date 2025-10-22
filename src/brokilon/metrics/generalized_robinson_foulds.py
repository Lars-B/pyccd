from brokilon.ccd.clades import DemeClade


def _get_deme_clades(tree, annotation_str):
    clades = set()
    for node in tree.traverse("levelorder"):
        # if len(node) > 1:  # this would ignore all leaf and internal nodes
        if len(node.children) != 1:  # this only ingores internal nodes, technically we care
            # about leaf nodes with this distance as they could contribute to it.
            cur_clade = {int(leaf.name) for leaf in node}
            clades.add(DemeClade(frozenset(cur_clade), deme=getattr(node, annotation_str)))
    return clades


def deme_robinson_foulds(tree_t1, tree_t2, annotation_str):
    # this is half as fast as the function above (which makes sense, its based on a weird function)
    clades1 = _get_deme_clades(tree_t1, annotation_str)
    clades2 = _get_deme_clades(tree_t2, annotation_str)
    return len(clades1.symmetric_difference(clades2))


def deme_robinson_foulds_reference(tree, reference_tree, annotation_str):
    # Essentially this is half the RF distance, all clades in the reference that are not in the tree
    clades_ref = _get_deme_clades(reference_tree, annotation_str)
    clades_tree = _get_deme_clades(tree, annotation_str)
    # All clades that are in the reference but not in the given tree (i.e. wrong clades if
    # reference is the truth)
    return len(clades_ref.difference(clades_tree))


if __name__ == '__main__':
    from pathlib import Path
    from brokilon.core import read_nexus_trees

    tree_file = (
        f"{Path(__file__).parent.absolute().parent.parent.parent}"
        f"/examples/data/h3n2-bdmm.h3n2_2deme.trees"
    )
    trees = read_nexus_trees(tree_file)
    t1 = trees[10]
    t2 = trees[11]
    d = deme_robinson_foulds(t1, t2, annotation_str="type")
    from brokilon.metrics import robinson_foulds, robinson_foulds_ccd

    regular_rf = robinson_foulds(t1, t2)
    ccd_rf = robinson_foulds_ccd(t1, t2)

    print(f"Annotated: {d}  ||| Regular: {regular_rf} ||| CCD: {ccd_rf}")
