from collections import defaultdict
from collections import deque
from statistics import mean

from brokilon.ccd.clades import DemeClade
from brokilon.ccd.types import CladeSplitInfo
from brokilon.core import Tree


def _find_deme(node, geo_ann_str):
    """
    From the node find the nearest node below that is either a leaf or a binary node.
    This is to ignore nodes that are along edges with state changes between binary nodes of the
    tree. (Might be possible to summarize the number of changes at some point but not now!)

    :param node: A node in a tree
    :param geo_ann_str: The attribute name to return
    :return: returns the value of the attribute of the nearest binary/leaf node
    """
    if node.is_leaf() or len(node.children) == 2:
        return getattr(node, geo_ann_str)

    # find the closest binary node
    cur_node = node.children[0]
    while len(cur_node.children) == 1:
        cur_node = cur_node.children[0]
    return getattr(cur_node, geo_ann_str)


def get_geo_map(trees, geo_ann_str, ccd_type=1):
    clade_count_map = defaultdict(int)
    clade_split_count_map = defaultdict(lambda: defaultdict(int))
    branch_lenghts_map = defaultdict(list)

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

            c0_deme = _find_deme(node.children[0], geo_ann_str)
            c1_deme = _find_deme(node.children[1], geo_ann_str)

            c0_clade = DemeClade(frozenset(c0_leaves), deme=c0_deme)
            c1_clade = DemeClade(frozenset(c1_leaves), deme=c1_deme)

            clade_count_map[parent_clade] += 1
            branch_lenghts_map[parent_clade].append(node.dist)
            if min(c0_leaves) < min(c1_leaves):
                clade_split_count_map[parent_clade][(c0_clade, c1_clade)] += 1
            else:
                clade_split_count_map[parent_clade][(c1_clade, c0_clade)] += 1
        elif node.is_leaf():
            leaf_clade = DemeClade(
                frozenset({int(node.name)}),
                deme=getattr(node, geo_ann_str)
            )
            branch_lenghts_map[leaf_clade].append(node.dist)

    if ccd_type == 1:
        # convert counts to probabilities
        for clade in clade_count_map:
            for split in clade_split_count_map[clade]:
                # todo this is where we could add a log version for precision and overflow.
                clade_split_count_map[clade][split] /= clade_count_map[clade]
        # convert to a dict to avoid defaultdict behaviour downstream
        return ({clade: dict(splits) for clade, splits in clade_split_count_map.items()},
                branch_lenghts_map,
                clade_count_map)
    elif ccd_type == 0:
        n_trees = len(trees)

        # todo missing expansion for now.. not well defined for these?

        geo_ccd0_probabilities = {}

        def recursive_prob_computer(clade):
            if clade in geo_ccd0_probabilities:
                return geo_ccd0_probabilities[clade]
            clade_value = clade_count_map.get(clade, 0) / n_trees

            # leaf
            if len(clade) == 1:
                geo_ccd0_probabilities[clade] = 1.0
                return 1.0

            phygeo_ccd0_map[clade] = {}
            part_products = []
            for left, right in clade_split_count_map[clade]:
                left_prob = recursive_prob_computer(left)
                right_prob = recursive_prob_computer(right)
                product = left_prob * right_prob
                part_products.append(product)
            total = sum(part_products)

            if total > 0:
                for (left, right), prod in zip(clade_split_count_map[clade], part_products):
                    phygeo_ccd0_map[clade][(left, right)] = prod / total
            else:
                raise ValueError("Need to implement a log fallback option")

            clade_prob = clade_value * total
            geo_ccd0_probabilities[clade] = clade_prob

            child_probs = defaultdict(float)
            queue = deque()

            visited = set()

            # assigning the frequency of each root to it as a starting probability
            for clade, clade_count in clade_count_map.items():
                if len(clade) == len(trees[0]):
                    child_probs[clade] = clade_count / n_trees
                    queue.append(clade)
                    visited.add(clade)

            while queue:
                clade = queue.popleft()
                parent_prob = child_probs[clade]

                for (left, right), ccp in phygeo_ccd0_map.get(clade, {}).items():
                    child_probs[left] += parent_prob * ccp
                    child_probs[right] += parent_prob * ccp

                    for child in (left, right):
                        if child in phygeo_ccd0_map and child not in visited:
                            queue.append(child)
                            visited.add(child)
            return clade_prob

        phygeo_ccd0_map = {}
        for clade in clade_count_map:
            recursive_prob_computer(clade)
        return phygeo_ccd0_map, branch_lenghts_map, clade_count_map
    else:
        raise ValueError(f"Unknown ccd_type: {ccd_type}")


def get_geo_map_tree(geo_ccd_map, geo_ann_str, taxon_map=None,
                     branch_length_map={}, clade_count_map=None):
    seen_resolved_clades = {}
    for clade in sorted(geo_ccd_map.keys(), key=len):
        if len(clade) == 2:
            best_split, prob = max(geo_ccd_map[clade].items(), key=lambda item: item[1])
            # todo handling ties somehow in the MAP tree
            #  (see transmission stuff)
            # ties = [k for k, v in geo_ccd_map[clade].items() if v == prob]
            # if len(ties):
            # print("TIEBREAKING FOR A LEAF IN EFFECT.")
            seen_resolved_clades[clade] = CladeSplitInfo(best_split, prob)
        else:
            for split, prob in geo_ccd_map[clade].items():
                split_c1, split_c2 = split
                if len(split_c1) == 1:
                    # This is fine if the leafs are all annotated with the same state,
                    # otherwise this needs to be changed to
                    # observed leaf with specific annotation / number of trees
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
    max_key_value = len(max(seen_resolved_clades.keys()))
    all_root_clades = [k for k in seen_resolved_clades if len(k) == max_key_value]
    max_root_prob = max(seen_resolved_clades[root].prob for root in all_root_clades)
    best_roots = [root for root in all_root_clades
                  if seen_resolved_clades[root].prob == max_root_prob]

    # TODO This would be a good situation to compute all MAP trees, for all roots?

    if len(best_roots) > 1:
        if clade_count_map:
            most_freq_root = max(best_roots, key=lambda x: clade_count_map[x])
            if not any(
                    [clade_count_map[root] == clade_count_map[most_freq_root]
                     for root in best_roots if root != most_freq_root]
            ):
                # checking if there are any other root clades with the same sample frequency
                working_list = [most_freq_root]
            else:
                raise ValueError("The root clade can not be resolved, there are root caldes"
                                 "with the same CCD probability and the same sampling frequency."
                                 "There is not method implemented to resolve such a conflict.")
        else:
            raise NotImplementedError("There is no way to resolve the root clade,"
                                      "multiple roots have the same CCD probability,"
                                      "you can pass the clade_count_map to resolve this based"
                                      "on the sampling frequency of a clade.")
    else:
        # only one root with  the highest probability, easy case
        working_list = [best_roots[0]]
    while working_list:
        cur_parent = working_list.pop()
        split, prob = seen_resolved_clades[cur_parent]
        output[cur_parent] = split
        working_list.extend(child for child in split if len(child) > 1)
    return get_tree_from_dict_of_splits(
        output, geo_ann_str, geo_ccd_map, taxon_map, branch_length_map
    )


def get_tree_from_dict_of_splits(
        splits, geo_ann_str,
        geo_ccd_map,
        taxon_map=None,
        branch_length_map={}
):
    root = max(splits.keys())

    output_tree = Tree(support=geo_ccd_map[root][splits[root]],
                       dist=mean(branch_length_map[root]) if branch_length_map else 0,
                       name="root")

    output_tree.add_feature(geo_ann_str, root.deme)
    icount = 1

    def recursive_children(node, new_split):
        nonlocal icount, splits, branch_length_map, geo_ccd_map
        c1, c2 = new_split

        def add_clade(node, clade):
            nonlocal icount, splits, geo_ccd_map
            if len(clade) == 1:
                label = taxon_map[next(iter(clade.clade))] if taxon_map else next(iter(clade.clade))
                leaf = node.add_child(
                    name=str(label),
                    dist=mean(branch_length_map.get(clade, [1])),
                    # Currently 1.0 because we assume leafs are always fixed in their annotation
                    # geo_ccd_map[clade][splits[clade]]
                    support=1.0
                )
                leaf.add_feature(geo_ann_str, clade.deme)
            else:
                internal_node = node.add_child(
                    name=f"internal_{icount}",
                    dist=mean(branch_length_map.get(clade, [1])),
                    support=geo_ccd_map[clade][splits[clade]]
                )
                internal_node.add_feature(geo_ann_str, clade.deme)
                icount += 1
                recursive_children(internal_node, splits[clade])

        add_clade(node, c1)
        add_clade(node, c2)

    recursive_children(output_tree, splits[root])

    return output_tree


def get_all_trees_represented(geo_ccd_map):
    # todo make a count only function too....?
    working_list_all_trees = {}

    for clade in sorted(geo_ccd_map.keys(), key=len):

        if len(clade) == 2:
            # assuming the leaf labels are fixed,
            # otherwise this probably will need to be handled differently!!!
            # working_list_all_trees[clade] = (
            #     [k for k in geo_ccd_map[clade].keys()]
            # )
            # todo work out the cherry case to build trees properly...
            #  best to add them here to avoid multiple loopings....
            continue
        elif len(clade) == 3:
            all_current_splits = [k for k in geo_ccd_map[clade].keys()]
            working_list_all_trees[clade] = (
                [split for split in all_current_splits]
            )
        else:
            all_current_splits = [k for k in geo_ccd_map[clade].keys()]

            for child1, child2 in all_current_splits:
                current_root = {clade: (child1, child2)}
                if len(child1) > 2:
                    all_c1_subtrees = working_list_all_trees[child1]
                else:
                    # todo this is wrong...
                    all_c1_subtrees = [k for k in geo_ccd_map[child2].keys()]

                partial_tree_list = []
                for c1split in all_c1_subtrees:
                    extended = current_root.copy()
                    extended[child1] = c1split
                    partial_tree_list.append(extended)

                if len(child2) > 2:
                    all_c2_subtrees = working_list_all_trees[child2]
                else:
                    # todo this is wrong...
                    all_c2_subtrees = [k for k in geo_ccd_map[child2].keys()]

                full_tree_list = []

                for c2split in all_c2_subtrees:
                    for partial in partial_tree_list:
                        to_be_made_full = partial.copy()
                        to_be_made_full[child2] = c2split
                        full_tree_list.append(to_be_made_full)
                if clade in working_list_all_trees:
                    working_list_all_trees[clade].extend(full_tree_list)
                else:
                    working_list_all_trees[clade] = full_tree_list

    n_taxa = max(len(k) for k in working_list_all_trees)
    all_roots = [k for k in working_list_all_trees if len(k) == n_taxa]
    # todo missing cherries but otherwise we should be able to call this:
    example_tree = (
        get_tree_from_dict_of_splits(working_list_all_trees[all_roots[1]][1])
    )

    print("help me...")
    # todo convert to trees and output list of trees.
    # todo implement function to just count number of trees,
    #  should be much easier...
    return []


if __name__ == '__main__':
    from pathlib import Path
    from brokilon.core.read_nexus import read_nexus_trees

    tree_file = (f"{Path(__file__).parent.absolute().parent.parent.parent.parent}/examples/"
                 f"data/h3n2-bdmm.h3n2_2deme.trees")
    trees, taxon_map = read_nexus_trees(tree_file, parse_taxon_map=True)

    trees = trees[int(len(trees) * 0.1):]
    # geo_ccd_map = get_geo_map(trees, geo_ann_str="type", ccd_type=1)
    geo_ccd_map, branch_length_map, _ = get_geo_map(trees, geo_ann_str="type", ccd_type=0)

    map_tree = get_geo_map_tree(geo_ccd_map, geo_ann_str="type",
                                taxon_map=taxon_map,
                                branch_length_map=branch_length_map)
    print(map_tree.write(format=5, features=["type"], format_root_node=True))
