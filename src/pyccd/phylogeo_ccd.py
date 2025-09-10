from collections import defaultdict
from dataclasses import dataclass

from pyccd.ccd0_attempt import CladeSplitInfo
from pyccd.transmission_ccd import BaseClade
from pyccd.tree import Tree


@dataclass(frozen=True)
class DemeClade(BaseClade):
    """
    Clade with a compartment annotation (deme) that indicates a location or similar

    Attributes:
        deme (string): Indicates a compartment this clade belongs to.
    """
    deme: str


def get_geo_map(trees, geo_ann_str, ccd_type=1):
    # todo unify to be ccd1 and ccd0... ccd_type WIP
    clade_count_map = defaultdict(int)
    clade_split_count_map = defaultdict(lambda: defaultdict(int))

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
            # non leaf nodes
            c0_leaves = {int(leaf.name) for leaf in node.children[0]}
            c1_leaves = {int(leaf.name) for leaf in node.children[1]}
            parent_leaves = frozenset(c0_leaves.union(c1_leaves))

            parent_clade = DemeClade(parent_leaves, deme=getattr(node, geo_ann_str))

            # todo start of cleanup
            # todo this needs a cleanup, maybe a for loop
            if node.children[0].is_leaf() or len(node.children[0].children) == 2:
                c0_deme = getattr(node.children[0], geo_ann_str)
            else:
                # find closest binary node
                cur_node = node.children[0].children[0]
                while len(cur_node.children) == 1:
                    cur_node = cur_node.children[0]
                c0_deme = getattr(cur_node, geo_ann_str)

            if node.children[1].is_leaf() or len(node.children[1].children) == 2:
                c1_deme = getattr(node.children[1], geo_ann_str)
            else:
                # find closest binary node...
                cur_node = node.children[1].children[0]
                while len(cur_node.children) == 1:
                    cur_node = cur_node.children[0]
                c1_deme = getattr(cur_node, geo_ann_str)
            # todo end of cleanup


            c0_clade = DemeClade(frozenset(c0_leaves), deme=c0_deme)
            c1_clade = DemeClade(frozenset(c1_leaves), deme=c1_deme)

            clade_count_map[parent_clade] += 1
            if min(c0_leaves) < min(c1_leaves):
                clade_split_count_map[parent_clade][(c0_clade, c1_clade)] += 1
            else:
                clade_split_count_map[parent_clade][(c1_clade, c0_clade)] += 1
        # todo do we need the leaves here? I don't think so but maybe later?...
        #  for the current way, the leaf annotations are set so all the same...
        # elif len(node) == 1:
        #    # leaf node
        #    geo_ann = getattr(node, geo_ann_str)
        #    clade = frozenset(int(node.name))
        #    leaf_clade = DemeClade(deme=geo_ann, clade=clade)
        #
        # else:
        #     raise ValueError("Something is wrong with this tree!")

    if ccd_type == 1:
        # convert counts to probabilities
        for clade in clade_count_map:
            for split in clade_split_count_map[clade]:
                # todo this is where we could add a log version for precision and overflow.
                clade_split_count_map[clade][split] /= clade_count_map[clade]
    elif ccd_type == 0:
        raise NotImplementedError("WIP")
    else:
        raise ValueError(f"Unknown ccd_type: {ccd_type}")

    # convert to a dict to avoid defaultdict behaviour downstream
    return {clade: dict(splits) for clade, splits in clade_split_count_map.items()}


def get_geo_map_tree(geo_ccd_map, geo_ann_str, taxon_map):
    seen_resolved_clades = {}
    for clade in sorted(geo_ccd_map.keys(), key=len):
        if len(clade) == 2:
            best_split, prob = max(geo_ccd_map[clade].items(), key=lambda item: item[1])
            ties = [k for k, v in geo_ccd_map[clade].items() if v == prob]
            if len(ties):
                print("TIEBREAKING FOR A LEAF IN EFFECT.")
            seen_resolved_clades[clade] = CladeSplitInfo(best_split, prob)
        else:
            for split, prob in geo_ccd_map[clade].items():
                split_c1, split_c2 = split
                if len(split_c1) == 1:
                    c1_prob = 1.0
                else:
                    c1_prob = seen_resolved_clades[split_c1].prob

                if len(split_c2) == 1:
                    c2_prob = 1.0
                else:
                    c2_prob = seen_resolved_clades[split_c2].prob
                this_split_probability = (c1_prob * c2_prob * prob)

                if (clade not in seen_resolved_clades or
                        seen_resolved_clades[clade].prob <= this_split_probability):
                    # todo think about tie breaking here, would be the prob == cases
                    seen_resolved_clades[clade] = CladeSplitInfo(split, prob)

    output = {}
    all_root_clades = [k for k in seen_resolved_clades.keys()
                       if len(k) == len(max(seen_resolved_clades.keys()))]
    max_root_prob = max(seen_resolved_clades[root].prob for root in all_root_clades)
    best_roots = [root for root in all_root_clades
                  if seen_resolved_clades[root].prob == max_root_prob]
    if len(best_roots) > 1:
        raise NotImplementedError("WIP: More than one root clade has max probability, "
                                  "not implemented...")

    working_list = [best_roots[0]]
    while working_list:
        cur_parent = working_list.pop()
        split, prob = seen_resolved_clades[cur_parent]
        output[cur_parent] = split
        working_list.extend(child for child in split if len(child) > 1)
    return get_tree_from_dict_of_splits(output, geo_ann_str, taxon_map)


def get_tree_from_dict_of_splits(splits, geo_ann_str, taxon_map):
    root = max(splits.keys())
    output_tree = Tree(support=0, dist=0, name="root")
    output_tree.add_feature(geo_ann_str, root.deme)
    icount = 1

    def recursive_children(node, new_split):
        nonlocal icount, splits
        c1, c2 = new_split

        def add_clade(node, clade):
            nonlocal icount, splits
            if len(clade) == 1:
                label = taxon_map[next(iter(clade.clade))]
                leaf = node.add_child(name=str(label), dist=1)
                leaf.add_feature(geo_ann_str, clade.deme)
            else:
                internal_node = node.add_child(name=f"internal_{icount}", dist=1)
                internal_node.add_feature(geo_ann_str, clade.deme)
                icount += 1
                recursive_children(internal_node, splits[clade])
        add_clade(node, c1)
        add_clade(node, c2)

    recursive_children(output_tree, splits[root])

    return output_tree


if __name__ == '__main__':
    from pathlib import Path
    from pyccd.read_nexus import read_nexus_trees

    tree_file = (f"{Path(__file__).parent.absolute().parent.parent}/examples/"
                 f"data/h3n2-bdmm.h3n2_2deme.trees")
    trees, taxon_map = read_nexus_trees(tree_file,
                                        breath_trees=False,
                                        label_transm_history=False,
                                        parse_taxon_map=True)

    # todo refactor and make this work for ccd0 and ccd1 to replace the old code
    # todo preliminary branch lengths would be nice for visualization stuff...

    trees = trees[int(len(trees) * 0.1):]
    geo_ccd_map = get_geo_map(trees, geo_ann_str="type")
    map_tree = get_geo_map_tree(geo_ccd_map, geo_ann_str="type", taxon_map=taxon_map)
    print(map_tree.write(format=5, features=["type"], format_root_node=True))
