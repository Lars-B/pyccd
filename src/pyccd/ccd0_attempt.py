from collections import defaultdict
from collections import namedtuple
from math import log
from typing import Any

from pyccd.tree import Tree

CladeSplitInfo = namedtuple("CladeSplitInfo", ["split", "prob"])


def expand(observed_clades, observed_clade_splits):
    expanded_clade_partitions = defaultdict(set)

    for clade, splits in observed_clade_splits.items():
        for left, right in splits:
            if min(left) < min(right):
                expanded_clade_partitions[clade].add((left, right))
            else:
                expanded_clade_partitions[clade].add((right, left))

    clade_buckets = defaultdict(set)
    for c in observed_clades:
        clade_buckets[len(c)].add(c)

    for clade in observed_clades:
        n = len(clade)
        if n <= 2:
            continue  # skip leaves and cherries

        existing_splits = set(expanded_clade_partitions.get(clade, set()))
        # clade_list = sorted(clade)

        # consider only left sizes from 1 to n//2
        for left_size in range(1, n // 2 + 1):
            for left in clade_buckets[left_size]:
                if left.issubset(clade):
                    right = clade - left
                    if right in clade_buckets[len(right)]:
                        # canonical ordering
                        l, r = (left, right) if min(left) < min(right) else (right, left)
                        split = (l, r)
                        if split not in existing_splits:
                            expanded_clade_partitions[clade].add(split)
                            # print(f"Found new split for {clade}: {split}")

    return expanded_clade_partitions


def get_maps_full(trees: list[Tree]) \
        -> tuple[defaultdict[Any, int], defaultdict[Any, int]]:
    """
    From a list of trees, return relevant CCD maps from clades/clade splits to counts.

    :param trees: list of input trees
    :return: maps for CCDs, clades to occurrences (m1), clades to clade splits (m2), unique trees
    """
    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled

    for ix, t in enumerate(trees):

        for node in t.traverse("levelorder"):
            parent_clade = frozenset(int(leaf.name) for leaf in node)
            m1[parent_clade] += 1
            if node.children:
                child0_clade = frozenset(sorted(int(leaf.name) for leaf in node.children[0]))
                child1_clade = frozenset(sorted(int(leaf.name) for leaf in node.children[1]))
                m2[(parent_clade, child0_clade, child1_clade)] += 1
    return m1, m2


def get_ccd0(trees):
    n_trees = len(trees)

    m1, m2 = get_maps_full(trees)

    clade_partitions = defaultdict(list)

    for (parent, child1, child2), count in m2.items():
        clade_partitions[parent].append((child1, child2))

    # todo what is happening here, put this together with get maps...
    expanded_test = expand(set(m1.keys()), clade_partitions)
    clade_partitions = expanded_test

    ccd0_probabilities = {}

    # todo need to extract this funciton...
    def recursive_prob_computer(clade):
        nonlocal ccd0_probabilities

        if clade in ccd0_probabilities:
            return ccd0_probabilities[clade]

        clade_value = m1.get(clade, 0) / n_trees

        # leaf
        if len(clade) == 1:
            ccd0_probabilities[clade] = 1.0
            return 1.0

        # cherry
        if len(clade) == 2:
            # [(k,v) for k,v in ccd0_probabilities.items() if len(k) == 2]
            assert len(clade_partitions[clade]) == 1, ("Cherry split in more than one way? "
                                                       "impossible!")
            left, right = list(clade_partitions[clade])[0]
            # recursion for the two leaves shouldn't do much?...
            recursive_prob_computer(left)
            recursive_prob_computer(right)
            ccd0_probabilities[clade] = clade_value
            return clade_value

        # general clade
        # todo fix this to be properly computed...
        partitions_ccp[clade] = {}
        part_products = []
        for left, right in clade_partitions[clade]:
            left_probability = recursive_prob_computer(left)
            right_probability = recursive_prob_computer(right)
            product = (left_probability * right_probability)
            part_products.append(product)

        total = sum(part_products)

        if total > 0:
            for (left, right), prod in zip(clade_partitions[clade], part_products):
                partitions_ccp[clade][(left, right)] = prod / total
        else:
            raise ValueError("Need to implement this fall back...")

        clade_prob = clade_value * total
        ccd0_probabilities[clade] = clade_prob
        # print(f"    Clade: {sorted(clade)} CCD0 Probability: {clade_prob}")

        child_probabilities = defaultdict(float)
        root_clade = max(m1.keys())
        child_probabilities[root_clade] = 1.0  # root probability starts at 1

        from collections import deque
        queue = deque([root_clade])

        while queue:
            clade = queue.popleft()
            parent_prob = child_probabilities[clade]

            for (left, right), ccp in partitions_ccp.get(clade, {}).items():
                # propagate probability to children
                child_probabilities[left] += parent_prob * ccp
                child_probabilities[right] += parent_prob * ccp

                # add children to the queue if they have further partitions
                if left in partitions_ccp:
                    queue.append(left)
                if right in partitions_ccp:
                    queue.append(right)

        return clade_prob

    # todo rename the variables to be more consistent across the board...
    # using the recursion form here on...
    partitions_ccp = {}
    for clade in m1:
        recursive_prob_computer(clade)

    return partitions_ccp


def get_tree_probability(tree, partitions_ccp, use_log=False):
    prob = 0.0 if use_log else 1
    for node in tree.traverse("levelorder"):
        if len(node) > 2:
            # Cherry or bigger
            left_clade = frozenset(sorted(int(leaf.name) for leaf in node.children[0]))
            right_clade = frozenset(sorted(int(leaf.name) for leaf in node.children[1]))
            split = (left_clade, right_clade) if min(left_clade) < min(right_clade) \
                else (right_clade, left_clade)
            parent_clade = left_clade | right_clade
            ccp = partitions_ccp.get(parent_clade, {}).get(split, 0.0)
            if ccp == 0.0:
                prob = float("-inf") if use_log else 0.0
                break
            if use_log:
                prob += log(ccp)
            else:
                prob *= ccp
    return float(prob)


def get_map_tree(partitions_ccp):
    all_clades_sorted = sorted(partitions_ccp.keys(), key=len)

    seen_resolved_clades = {}
    for cur_clade in all_clades_sorted:
        assert len(cur_clade) > 2, "Cherries are not relevant for this?"
        if len(cur_clade) <= 3:
            # for a triplet we need to simply choose the max
            cur_clade_split, cur_clade_split_prob = max(partitions_ccp[cur_clade].items(),
                                                        key=lambda item: item[1])
            assert cur_clade not in seen_resolved_clades, "Triplets should only get into here once"
            seen_resolved_clades[cur_clade] = CladeSplitInfo(cur_clade_split, cur_clade_split_prob)
        else:
            for cur_clade_split, cur_clade_split_prob in partitions_ccp[cur_clade].items():
                cur_clade_child1, cur_clade_child2 = cur_clade_split

                if len(cur_clade_child1) < 3:
                    cur_clade_child1_prob = 1.0
                else:
                    assert cur_clade_child1 in seen_resolved_clades, "This should be impossible"
                    cur_clade_child1_prob = seen_resolved_clades[cur_clade_child1].prob

                if len(cur_clade_child2) < 3:
                    cur_clade_child2_prob = 1.0
                else:
                    assert cur_clade_child2 in seen_resolved_clades, "This should be impossible"
                    cur_clade_child2_prob = seen_resolved_clades[cur_clade_child2].prob

                cur_clade_split_map_probabiltiy = cur_clade_child1_prob * cur_clade_child2_prob * cur_clade_split_prob
                if cur_clade in seen_resolved_clades:
                    if seen_resolved_clades[cur_clade].prob <= cur_clade_split_map_probabiltiy:
                        seen_resolved_clades[cur_clade] = CladeSplitInfo(cur_clade_split,
                                                                         cur_clade_split_map_probabiltiy)
                else:
                    seen_resolved_clades[cur_clade] = CladeSplitInfo(cur_clade_split,
                                                                     cur_clade_split_map_probabiltiy)

    output = {}
    working_list = [max(seen_resolved_clades.keys())]

    while working_list:
        cur_parent = working_list.pop()
        split, prob = seen_resolved_clades[cur_parent]
        output[cur_parent] = split
        for split_child in split:
            if len(split_child) > 2:
                working_list.append(split_child)

    return output


def get_tree_from_dict_of_splits(splits):
    # todo need to fix the distances/branch lengths to make an ultrametric tree
    output_tree = Tree(support=0, dist=0, name="root")
    icount = 1

    def recursive_children(node, new_split):
        nonlocal splits, icount
        clade_c1, clade_c2 = new_split

        if len(clade_c1) == 2:
            # add two leaf nodes
            l1, l2 = clade_c1
            internal = node.add_child(name=f"internal_{icount}", dist=1)
            icount += 1
            internal.add_child(name=f"leaf_{l1}", dist=1)
            internal.add_child(name=f"leaf_{l2}", dist=1)
        elif len(clade_c1) == 1:
            # add a single leaf node
            label = next(iter(clade_c1))
            node.add_child(name=f"leaf_{label}", dist=1)
        else:
            c1 = node.add_child(name=f"child_internal_{icount}", dist=1)
            icount += 1
            recursive_children(c1, splits[clade_c1])

        if len(clade_c2) == 2:
            # add two leaf nodes
            l1, l2 = clade_c2
            internal = node.add_child(name=f"internal_{icount}", dist=1)
            icount += 1
            internal.add_child(name=f"leaf_{l1}", dist=1)
            internal.add_child(name=f"leaf_{l2}", dist=1)
        elif len(clade_c2) == 1:
            # add a single leaf node
            label = next(iter(clade_c2))
            node.add_child(name=f"leaf_{label}", dist=1)
        else:
            # recursion further down
            c2 = node.add_child(name=f"child_internal_{icount}", dist=1)
            icount += 1
            recursive_children(c2, splits[clade_c2])

    recursive_children(output_tree, splits[max(splits.keys())])
    print(output_tree.write(format=5))
    return None


if __name__ == '__main__':
    from pathlib import Path

    java_tree_probs = {}

    from pyccd.read_nexus import read_nexus_trees

    tree_file = f"{Path(__file__).parent.absolute().parent.parent}/examples/data/30Taxa.trees"
    trees, taxon_map = read_nexus_trees(tree_file,
                                        breath_trees=False,
                                        label_transm_history=False,
                                        parse_taxon_map=True)

    trees = trees[int(len(trees) * 0.1):]  # burn in deletion....
    partitions_ccp = get_ccd0(trees)

    output = get_map_tree(partitions_ccp)

    get_tree_from_dict_of_splits(output)

    # todo compare ccd0 MAP to the treeannotator output
