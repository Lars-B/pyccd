from brokilon.ccd.clades import DemeClade


def _get_deme_clades(tree, annotation_str):
    clades = set()
    for node in tree.traverse("levelorder"):
        if len(node) > 1:
            cur_clade = {int(leaf.name) for leaf in node}
            clades.add(DemeClade(frozenset(cur_clade), deme=getattr(node, annotation_str)))
    return clades


def deme_robinson_foulds(tree_t1, tree_t2, annotation_str):
    # this is half as fast as the function above (which makes sense, its based on a weird function)
    clades1 = _get_deme_clades(tree_t1, annotation_str)
    clades2 = _get_deme_clades(tree_t2, annotation_str)
    return len(clades1.symmetric_difference(clades2))


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
    from brokilon.metrics import robinson_foulds

    regular_rf = robinson_foulds(t1, t2)

    print(f"Annotated: {d}  ||| Regular: {regular_rf}")
