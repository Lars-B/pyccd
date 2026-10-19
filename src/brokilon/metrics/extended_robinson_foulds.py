from brokilon.ccd.clades import DemeClade
from brokilon.ccd.domain.transmission.transmission_ccd import _sanitize_transm_ancest


def _get_deme_clades(tree, annotation_str):
    clades = set()
    for node in tree.traverse("levelorder"):
        # if len(node) > 1:  # this would ignore all leaf and internal nodes
        if len(node.children) != 1:  # this only ignores internal nodes, technically we care
            # about leaf nodes with this distance as they could contribute to it.
            if not hasattr(node, annotation_str):
                try:
                    probs = [x for x in
                             getattr(node, f"{annotation_str}_set_prob").strip("[]\"\'").split(
                                 ",")]
                    types = [x for x in
                             getattr(node, f"{annotation_str}_set").strip("[]\"\'").split(",")]
                except Exception:
                    raise ValueError(f"Something is wrong with the tree file, likely inconsistent "
                                     f"labels from treeannotator which is infuriatingly stupid! "
                                     f"Double check whether your tree contains the"
                                     f" {annotation_str!r} attributes at nodes.")

                cur_clade = {int(leaf.name) for leaf in node}
                # this is a bypass, if there is no attribute of the required annotation it will
                # try to use the set and set prob to extract the most likely one, this is stupid
                # and all because treeannotoator doesn't annotate at least the leaves
                clades.add(DemeClade(frozenset(cur_clade), deme=types[probs.index(max(probs))]))
            else:
                cur_clade = {int(leaf.name) for leaf in node}
                clades.add(
                    DemeClade(
                        frozenset(cur_clade),
                        deme=getattr(node, annotation_str, "Missing")
                    )
                )
    return clades


def _get_transmission_clades(tree, annotation_str):
    clades = set()
    for node in tree.traverse("levelorder"):
        cur_clade = {int(leaf.name) for leaf in node}

        clades.add(DemeClade(
            frozenset(cur_clade),
            deme=_sanitize_transm_ancest(getattr(node, annotation_str, "Missing"))
        ))

    if all(clade.deme == "Missing" for clade in clades):
        raise ValueError(
            f"No nodes in the tree have the annotation {annotation_str!r}"
        )

    return clades


def _get_sranges_clades(tree, annotation_str):
    raise NotImplementedError


def extended_robinson_foulds(
        tree_t1, tree_t2, kind="deme", annotation_str=None
):
    # Clade extractor methods and their default annotation string
    extractors = {
        "deme": (_get_deme_clades, "type"),
        "transmission": (_get_transmission_clades, "transm_ancest"),
        "sranges": (_get_sranges_clades, None),
    }

    try:
        extractor, default_annotation = extractors[kind]
    except KeyError:
        raise ValueError(f"Unknown RF distance kind: {kind!r}") from None

    if annotation_str is None:
        annotation_str = default_annotation

    clades1 = extractor(tree_t1, annotation_str=annotation_str)
    clades2 = extractor(tree_t2, annotation_str=annotation_str)

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
    d = extended_robinson_foulds(t1, t2, kind="deme")
    from brokilon.metrics import robinson_foulds, robinson_foulds_ccd

    regular_rf = robinson_foulds(t1, t2)
    ccd_rf = robinson_foulds_ccd(t1, t2)

    print(f"Annotated: {d}  ||| Regular: {regular_rf} ||| CCD: {ccd_rf}")
