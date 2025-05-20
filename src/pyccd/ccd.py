"""
Module for the Conditional Clade Distribution implementation using maps for clades and clade splits.
"""
import re
from collections import defaultdict
from decimal import Decimal
import ete3
from numpy import random, log


def get_clades(tree: ete3.Tree) -> set[frozenset[str]]:
    """
    Get all clades of a given tree

    :param tree: an ete3 input tree
    :return: set of clades
    """
    treestr = tree.write(format=9)
    clades = set()
    if not (treestr[0] == '(' and treestr[-2] == ')' and treestr[-1] == ';'):
        raise ValueError("Invalid tree string given! (no ';' at the end)")
    opend = []
    re_brackets = re.compile(r"\(|\)")
    for i in range(1, len(treestr) - 2):
        if treestr[i] == '(':
            opend.append(i)
        elif treestr[i] == ')':
            if not opend:
                raise ValueError("Invalid tree string given! (to many ')')")
            cur = treestr[opend[-1]:i]
            clades.add(frozenset(re.sub(re_brackets, '', cur).split(',')))
            del opend[-1]
    if opend:
        raise ValueError("Invalid tree string given! (to many '(')")
    return clades


def get_maps(trees: list[ete3.Tree]) \
        -> tuple[defaultdict[str, int], defaultdict[str, int], dict[int, list]]:
    """
    From a list of trees, return relevant CCD maps from clades/clade splits to counts.

    :param trees: list of ete3 input trees
    :return: maps for CCDs, clades to occurrences (m1), clades to clade splits (m2), unique trees
    """
    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled
    uniques = {}

    seen = {}

    for ix, t in enumerate(trees):
        if not frozenset(sorted(get_clades(t))) in seen:
            seen[frozenset(sorted(get_clades(t)))] = ix
            uniques[ix] = []
        else:
            uniques[seen[frozenset(sorted(get_clades(t)))]].append(ix)

        for node in t.traverse("levelorder"):
            if len(node) > 2:
                c = node.children
                c0_leafs = set()
                for leaf in c[0]:
                    c0_leafs.add(int(leaf.name))
                c1_leafs = set()
                for leaf in c[1]:
                    c1_leafs.add(int(leaf.name))
                parent_clade = frozenset(sorted(c0_leafs.union(c1_leafs)))
                m1[parent_clade] += 1
                if min(c0_leafs) < min(c1_leafs):
                    m2[(parent_clade, frozenset(c0_leafs))] += 1
                else:
                    m2[(parent_clade, frozenset(c1_leafs))] += 1
    return m1, m2, uniques


def get_tree_probability(tree, m1, m2, use_log=False):
    """
    Calculate the probability of a tree given the occurrences of its clade and clade splits.

    :param tree: input tree
    :param m1: CCD map for clades
    :param m2: CCD map for clade splits
    :param use_log: Whether to use log transform for probabilities
    :return: Probability of a tree
    """
    # getcontext().prec = 20
    probability = 0 if use_log else 1
    for node in tree.traverse("levelorder"):
        if len(node) > 2:
            c = node.children
            c0_leafs = set()
            for leaf in c[0]:
                c0_leafs.add(int(leaf.name))
            c1_leafs = set()
            for leaf in c[1]:
                c1_leafs.add(int(leaf.name))
            parent_clade = frozenset(sorted(c0_leafs.union(c1_leafs)))
            if m1[parent_clade] != 0:
                leaf_set = frozenset(c0_leafs) if min(c0_leafs) < min(c1_leafs) \
                    else frozenset(c1_leafs)
                m2_value = m2[(parent_clade, leaf_set)]
                m1_value = m1[parent_clade]
                if use_log:
                    probability += log(m2_value / m1_value)
                else:
                    probability *= Decimal(m2_value) / Decimal(m1_value)
    return float(probability)


def get_ccd_tree_bottom_up(m1, m2):
    """
    From the maps of clade counts and clade split counts perform a dynamic program to calculate
    the CCD MAP tree.

    :param m1: Map for clade counts
    :param m2: Map for clade split counts
    :return: the CCD MAP tree
    """
    # initialize with root clade, empty tree, probability 1
    # working_list = [([max(m1.keys())], [], 1)]

    seen_resolved_clades = {}

    all_clades = sorted(list(m1.keys()), key=len)
    for current_clade in all_clades:
        for split in [i for i in m2 if i[0] == current_clade]:
            # this for loop needs to find the best split of the current parent, if any exists!
            child1 = split[1]
            child2 = split[0].difference(split[1])

            c1_prob, c2_prob = 0, 0

            if len(child1) < 3:
                # length 2 or 1 gives probability 1
                c1_prob = 1
            else:
                if child1 in seen_resolved_clades:
                    c1_prob, _ = seen_resolved_clades[child1]
            if len(child2) < 3:
                # length 2 or 1 gives probability 1
                c2_prob = 1
            else:
                if child2 in seen_resolved_clades:
                    c2_prob, _ = seen_resolved_clades[child2]

            # if (c1_prob * c2_prob) < prob:
            #     # The current split will not lead to an improvement
            #     continue

            # cur_prob = m2[split] / m1[split[0]]  # Prob of current parent, given split
            # best probability of current parent with split
            # split_prob = c1_prob * c2_prob * cur_prob
            split_prob = c1_prob * c2_prob * (m2[split] / m1[split[0]])
            # math.isclose instead of == for float comparison 0.1+0.2 != 0.3
            # if math.isclose(split_prob, prob) or split_prob > prob:
            if split[0] in seen_resolved_clades:
                # parent was already found, do we need to update?
                if seen_resolved_clades[split[0]][0] <= split_prob:
                    # this split has better probability,
                    # therefore update the seen_resolved_clades with the better split of split[0]
                    seen_resolved_clades[split[0]] = (split_prob, child1)
            else:
                # we have not seen the parent before
                seen_resolved_clades[split[0]] = (split_prob, child1)

    output = []
    # root = max(seen_resolved_clades.keys())
    # working_list = [root]
    working_list = [max(seen_resolved_clades.keys())]

    while working_list:
        cur_parent = working_list.pop()
        output.append((cur_parent, seen_resolved_clades[cur_parent][1]))
        if len(seen_resolved_clades[cur_parent][1]) > 2:
            working_list.append(seen_resolved_clades[cur_parent][1])
        if len(cur_parent.difference(seen_resolved_clades[cur_parent][1])) > 2:
            working_list.append(cur_parent.difference(seen_resolved_clades[cur_parent][1]))

    return get_tree_from_list_of_splits(output)


def get_tree_from_list_of_splits(splits) -> str:
    """
    From a list of splits create the corresponding tree as a newick string

    :param splits: list of splits
    :return: newick string of tree
    """
    # this is dependent on the current structure of how the greedy list of splits is created!
    n_taxa = len(splits[0][0])
    dist = 1
    # support is the rank of the node...
    cur_t = ete3.Tree(support=0, dist=0, name=",".join([str(i) for i in
    sorted(splits[0][0])]))
    for parent, child1 in splits:
        node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(parent)]))[0]
        child2 = parent.difference(child1)
        if len(child1) == 1:
            node.add_child(name=",".join([str(i) for i in sorted(child1)]), support=n_taxa - 1)
        else:
            node.add_child(name=",".join([str(i) for i in sorted(child1)]), support=dist)
            dist += 1
        if len(child2) == 1:
            node.add_child(name=",".join([str(i) for i in sorted(child2)]), support=n_taxa - 1)
        else:
            node.add_child(name=",".join([str(i) for i in sorted(child2)]), support=dist)
            dist += 1
        # IF the children nodes have 2 taxa the correspoinding leafs need to be added here
        if len(child1) == 2:
            node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(child1)]))[0]
            node.add_child(name=list(sorted(child1))[0], support=n_taxa - 1)
            node.add_child(name=list(sorted(child1))[1], support=n_taxa - 1)
        if len(child2) == 2:
            node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(child2)]))[0]
            node.add_child(name=list(sorted(child2))[0], support=n_taxa - 1)
            node.add_child(name=list(sorted(child2))[1], support=n_taxa - 1)

    # setting the node distances to ranks, no meaning just to have a ranked tree
    for node in cur_t.iter_descendants("postorder"):
        node.dist = node.support - node.up.support

    return cur_t.write(format=5)


# function to add a single tree to the maps m1 and m2
# def add_centroid(centroid, m1, m2):
#     for node in centroid.tree.traverse("levelorder"):
#         if len(node) > 2:
#             c = node.children
#             c0_leafs = set()
#             for leaf in c[0]:
#                 c0_leafs.add(int(leaf.name))
#             c1_leafs = set()
#             for leaf in c[1]:
#                 c1_leafs.add(int(leaf.name))
#             parent_clade = frozenset(sorted(c0_leafs.union(c1_leafs)))
#             if m1[parent_clade] == 0:
#                 m1[parent_clade] += 1
#                 if min(c0_leafs) < min(c1_leafs):
#                     m2[(parent_clade, frozenset(c0_leafs))] += 1
#                 else:
#                     m2[(parent_clade, frozenset(c1_leafs))] += 1
#             else:
#                 # x = int(m1[parent_clade]/9)
#                 # m1[parent_clade] += x
#                 m1[parent_clade] += 1
#                 if min(c0_leafs) < min(c1_leafs):
#                     # m2[(parent_clade, frozenset(c0_leafs))] += x
#                     m2[(parent_clade, frozenset(c0_leafs))] += 1
#                 else:
#                     # m2[(parent_clade, frozenset(c1_leafs))] += x
#                     m2[(parent_clade, frozenset(c1_leafs))] += 1
#     return m1, m2


def sample_tree_from_ccd(m1, m2, n=1) -> list[ete3.Tree]:
    """
    Given a CCD with m1 and m2, this function samples n trees proportional to
    their probabilities from this CCD.

    :param m1: Count of clades
    :param m2: Count of clade splits
    :param n: number of trees to sample
    :return: List of sampled trees
    """
    # sample n trees from the CCD distribution, relative to its clade probabilities in each step
    samples = []

    for _ in range(n):
        cur_sample = []
        # cur_sample.append(max(m1))

        working_list = [max(m1)]
        while working_list:
            cur_clade = working_list.pop()
            possible_splits = [(list(i), m2[i]) for i in m2 if list(i)[0] == cur_clade]

            cur_sum = m1[cur_clade]  # same as sum([i[1] for i in next_splits])
            cur_p = [i[1] / cur_sum for i in possible_splits]

            chosen_split = random.choice([i[0][1] for i in possible_splits], p=cur_p)
            remainder_split = cur_clade.difference(chosen_split)
            if len(chosen_split) > 2:
                working_list.append((chosen_split))
            if len(remainder_split) > 2:
                working_list.append(remainder_split)
            cur_sample.append((cur_clade, chosen_split))
        samples.append(get_tree_from_list_of_splits(cur_sample))
    return samples


# def sample_logprob_from_ccd(m1, m2, n=1):
#     # NOTE: this is fairly inefficient for some reason,
#     #  may need change in the future, same for the other sampling function
#     # sample n trees from the CCD distribution, relative to its clade probabilities in each step
#     # samples = []
#     probabilities = []
#
#     for _ in range(n):
#         # cur_sample = []
#         cur_prob = 0
#         # cur_sample.append(max(m1))
#
#         working_list = [max(m1)]
#         while working_list:
#             cur_clade = working_list.pop()
#             possible_splits = [(list(i), m2[i]) for i in m2 if list(i)[0] == cur_clade]
#
#             cur_sum = m1[cur_clade]  # same as sum([i[1] for i in next_splits])
#             cur_p = [i[1] / cur_sum for i in possible_splits]
#
#             chosen_split = random.choice([i[0][1] for i in possible_splits], p=cur_p)
#             remainder_split = cur_clade.difference(chosen_split)
#             if len(chosen_split) > 2:
#                 working_list.append((chosen_split))
#             if len(remainder_split) > 2:
#                 working_list.append(remainder_split)
#             # cur_sample.append((cur_clade, chosen_split))
#             cur_prob += log(m2[(cur_clade, chosen_split)] / m1[cur_clade])
#
#         # samples.append(get_tree_from_list_of_splits(cur_sample))
#         probabilities.append(cur_prob)
#     return probabilities


def calc_entropy(m1, m2):
    """
    For a given CCD via the maps m1 and m2, this calculates the entropy for it using the
    fromular from Lewis et al.

    :param m1: Count of clades
    :param m2: Count of clade splits
    :return: Entropy as a float
    """
    h_dict = defaultdict(lambda: 0)
    for c in sorted(m1.keys(), reverse=False, key=len):
        # iterate over all clades, from small to large
        h_dict[c] = 0
        c_children = [k for k in m2.keys() if k[0] == c]
        for _, child in c_children:
            p = m2[(c, child)] / m1[c]
            # if len(c) == 3:
            #     # if c has only 3 taxa we don't need to go for childrens entropy
            #     h_dict[c] -= p * np.log(p)
            # else:
            # if c has more than 3 taxa use formula
            h_dict[c] -= p * (log(p) - h_dict[child] - h_dict[c.difference(child)])
    return h_dict[max(m1.keys())]
