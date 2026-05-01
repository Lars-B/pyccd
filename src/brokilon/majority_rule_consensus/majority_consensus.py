from collections import defaultdict

from brokilon.ccd.domain.phylogeography import DemeClade
from brokilon.core import Tree


def regular_mrc(trees):
    print(len(trees))
    raise NotImplementedError("WIP")


# def majority_consensus_annoated(trees):
#
#     m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees,
#                                                                        type_str="Ancestry")
#     major_thr = 0.5
#     clades_majority = [k for k, v in m1.items() if v / len(trees) > major_thr]
#     if len(clades_majority) < len(trees[0]) / 4:
#         print(f"No clades found with threshold: {major_thr}, will divide by 2 and try again")
#         major_thr = major_thr / 2
#         clades_majority = [k for k, v in m1.items() if v / len(trees) > major_thr]
#
#     # todo  m1 does not contain leaf clades, but they are now important so we need to keep that in
#     #  mind
#
#     print(len(clades_majority))
#     # print(clades_majority)
#     for c in clades_majority:
#         print(c)
#     return consensus


def phygeo_mrc(trees, geo_ann_str, major_thr=0.5):
    """
    Computes MRC tree for annotated trees (type is input value).
    Assumes leaves are fixed states.

    :param trees: list of trees
    :param geo_ann_str: annotation that should be used for clades
    :return: MRC tree with annotations
    """
    if not 0.0 < major_thr < 1.0:
        raise ValueError(f"Major threshold must be between 0 and 1, but got {major_thr}")

    # todo with very low values this might not be able to construct at tree... catch that.

    clade_count_map = defaultdict(int)
    for node in (node for t in trees for node in t.traverse("levelorder")):
        # if geo_ann_str not in node.features:
        if len(node.children) == 1:
            # the structured approaches have internal nodes, skip them for now, might be
            # interesting to think about them later...
            continue
        if not hasattr(node, geo_ann_str):
            raise ValueError(f"Node does not have a geographic annotation of name '{geo_ann_str}'."
                             f"Available are {node.features}, maybe change the input value...")

        if len(node) > 1:
            parent_leaves = frozenset(int(leaf.name) for leaf in node)
            parent_clade = DemeClade(parent_leaves, deme=getattr(node, geo_ann_str))
            clade_count_map[parent_clade] += 1

    # clades_majority = [k for k, v in clade_count_map.items() if v / len(trees) > major_thr]
    clades_majority = []
    clade_support = {}
    for k, v in clade_count_map.items():
        if v / len(trees) > major_thr:
            clades_majority.append(k)
            clade_support[k] = v / len(trees)

    working_subtrees = {}
    for l in trees[0]:
        t = Tree(support=1.0, dist=1.0, name=l.name)
        t.add_feature(geo_ann_str, getattr(l, geo_ann_str))
        working_subtrees[l.name] = t

    for i, c in enumerate(sorted(clades_majority, key=len)):
        new_node = Tree(
            name=f"clade_{i}",
            support=clade_support[c],
            dist=1.0,
        )
        new_node.add_feature(geo_ann_str, c.deme)

        for taxon in c.clade:
            skip = False
            node = working_subtrees.pop(str(taxon))

            if isinstance(node, str) and node.startswith("clade_"):
                if node in working_subtrees:
                    if (
                            isinstance(working_subtrees[node], str)
                            and
                            working_subtrees[node].startswith("clade_")
                    ):
                        working_subtrees.pop(node)
                        skip = True
                    else:
                        node = working_subtrees.pop(node)
                else:
                    skip = True
            if not skip:
                new_node.add_child(node)
                working_subtrees[node.name] = f"clade_{i}"
            working_subtrees[str(taxon)] = f"clade_{i}"

        working_subtrees[f"clade_{i}"] = new_node

    try:
        output = working_subtrees[f"clade_{len(clades_majority) - 1}"]
    except KeyError:
        raise ValueError("Something went wrong when constructing the MRC tree.")
    return output


def transmission_mrc(trees):
    print(len(trees))
    raise NotImplementedError("WIP")


if __name__ == '__main__':
    from pathlib import Path
    from brokilon.core.read_nexus import read_nexus_trees

    # trees_file = (f"{Path(__file__).parent.absolute().parent.parent}"
    #              f"/examples/data/breath32sim.trees")
    tree_file = (f"{Path(__file__).parent.absolute().parent.parent}"
                 f"/examples/data/"
                 # f"roetzer40.trees"
                 f"h3n2-bdmm.h3n2_2deme.trees")
    # todo burnin missing...

    trees, taxon_map = read_nexus_trees(
        tree_file,
        breath_trees=False,
        label_transm_history=False,
        parse_taxon_map=True
    )
    phygeo_mrc(trees, geo_ann_str="type")
    # majority_consensus_annoated(trees)
