from collections import defaultdict
from math import log
from typing import Any

from pyccd.tree import Tree


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


# todo extract this into the example part and automate properly...
def compare_to_java_tree_probs(java_tree_probs, partitions_ccp):
    python_probs = {}
    for i, t in enumerate(trees):
        python_probs[i] = get_tree_probability(t, partitions_ccp)

    tolerance = 1e-12

    print("Look for mismatches between java and python....")
    prob_problem = []
    for k in python_probs.keys():
        v1 = python_probs[k]
        v2 = java_tree_probs[k]
        if abs(v1 - v2) > tolerance:
            prob_problem.append((k, v1, v2))
    print(prob_problem)


if __name__ == '__main__':
    from pathlib import Path

    java_tree_probs = {}
    # todo this should be automated properly, maybe write a little beast app for this...
    #  or does this already exist if we use treestat? or something???
    with (open(f"{Path(__file__).parent.absolute().parent.parent}/examples/data/java_probs.out")
          as f):
        for line in f:
            if not line.startswith("Tree "):
                continue
            cur_index = int(line.split("=", 1)[0].split(" ")[1])
            cur_prob = float(line.split("=", 1)[1].strip())
            java_tree_probs[cur_index] = cur_prob

    from pyccd.read_nexus import read_nexus_trees

    tree_file = f"{Path(__file__).parent.absolute().parent.parent}/examples/data/30Taxa.trees"
    trees = read_nexus_trees(tree_file, breath_trees=False, label_transm_history=False)
    # trees = trees[int(len(trees)*0.1):]  # burn in deletion....
    partitions_ccp = get_ccd0(trees)

    compare_to_java_tree_probs(java_tree_probs, partitions_ccp)

    # todo next steps are to get a ccd0 map tree
    #  then compare that to the treeannotator output
    # todo think about how to do this for ccd with transmission histories?...
    #  then apply to simulation datasets etc...
    # todo evaluate the breath stuff for increasing amount of trees
    #  compare posterior wiw vs the summary trees and compare the scorse?
    # todo write a cli interface to calculate the ccd1 and ccd0 map tree with this
    # todo provide a example script that checks it creates the same as treeannotator
    # todo entropy calculation needs to be adopted... can it be general for all types?
