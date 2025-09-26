from collections import defaultdict

from brokilon.phylogeography import DemeClade
from brokilon.transmission_ccd import get_transmission_maps


def majority_consensus_annoated(trees):
    consensus = 0

    # todo figure out if this works...
    # here we only need a map of clades and how often that specific clade was observed
    # that should be a separate function as ccd maps should always be probailities? ...

    print(len(trees))
    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees,
                                                                       type_str="Ancestry")
    major_thr = 0.5
    clades_majority = [k for k, v in m1.items() if v / len(trees) > major_thr]
    if len(clades_majority) < len(trees[0]) / 4:
        print(f"No clades found with threshold: {major_thr}, will divide by 2 and try again")
        major_thr = major_thr / 2
        clades_majority = [k for k, v in m1.items() if v / len(trees) > major_thr]

    # todo expand the Maps m1, m2 to include leaf clades etc? or are they just ignored in M1 but
    #  present, similar to what I did for CCD0 attempt
    #  in M2?...

    # todo  m1 does not contain leaf clades, but they are now important so we need to keep that in
    #  mind

    print(len(clades_majority))
    # print(clades_majority)
    for c in clades_majority:
        print(c)
    return consensus


def phygeo_mrc(trees, geo_ann_str):
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

    major_thr = 0.5
    clades_majority = [k for k, v in clade_count_map.items() if v / len(trees) > major_thr]
    # todo for this we want to include the leaf clades too

    # todo build a tree from these clades...
    # 1. start with largest clade (make sure there is just one but should be fine)
    # Maybe start from the leafs and then go up, merge all cherries, merge all other things etc...
    # this isn't as trivial as I thought it would be...

    return None


if __name__ == '__main__':
    from pathlib import Path
    from brokilon.core.read_nexus import read_nexus_trees

    # tree_file = (f"{Path(__file__).parent.absolute().parent.parent}"
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
