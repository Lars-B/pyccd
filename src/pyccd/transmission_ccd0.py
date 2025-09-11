from collections import defaultdict

from pyccd.ccd0_attempt import CladeSplitInfo
from pyccd.transmission_ccd import get_transmission_maps
from pyccd.tree import Tree


def compatible(split, parent):
    l, r = split

    pancest = False if parent.transm_ancest.startswith("Unknown") else int(parent.transm_ancest)
    lancest = False if l.transm_ancest.startswith("Unknown") else int(l.transm_ancest)
    rancest = False if r.transm_ancest.startswith("Unknown") else int(r.transm_ancest)

    if lancest in r.clade and rancest in l.clade:
        # this doesn't work
        return False

    if not any([pancest, lancest, rancest]):
        # all unknown so this is fine...
        return True

    if pancest in l.clade:
        if lancest == pancest:
            if rancest == pancest or not rancest:
                return True
        return False
    if pancest in r.clade:
        if rancest == pancest:
            if lancest == pancest or not lancest:
                return True
        return False

    if not any([lancest, rancest]):
        # both left and right are unknown which is fine
        return True

    # here the parent ancestor is not in either left or right
    # left and right have to be compatible though...
    if lancest in l.clade:
        if not rancest or rancest == lancest:
            return True
        return False
    if rancest in r.clade:
        if not lancest or lancest == rancest:
            return True
        return False

    if rancest not in l.clade:
        return False
    if lancest not in r.clade:
        return False

    print(
        f"-----------------\n"
        f"No case applied, will be marked as not compatible...\n"
        f"{parent}\n"
        f"{split}\n"
        f"-----------------"
    )

    return False


def expand_tccd0(clade_partitions):
    expanded_clade_partitions = defaultdict(set)

    for clade, splits in clade_partitions.items():
        for left, right in splits:
            if min(left.clade) < min(right.clade):
                expanded_clade_partitions[clade].add((left, right))
            else:
                expanded_clade_partitions[clade].add((right, left))

    clade_buckets = defaultdict(set)
    for c in clade_partitions.keys():
        clade_buckets[len(c)].add(c)

    total = len(clade_partitions)

    for i, clade in enumerate(clade_partitions.keys(), 1):
        n = len(clade)

        # progress prints
        if i % 100 == 0 or i == total:
            print(f"[{i}/{total}] Working on clade of size {n}")

        if n < 2:
            # we now care about cherries, might even have to do leafs?...
            continue

        existing_splits = set(expanded_clade_partitions.get(clade, set()))

        for left_size in range(1, n // 2 + 1):
            # print(f"Working on {left_size}")
            # this if is early stopping if there is nothing to combine with
            if clade_buckets[n - left_size]:
                for left in clade_buckets[left_size]:
                    if left.clade.issubset(clade.clade):
                        rightCladeSet = clade.clade - left.clade
                        matches = [tc for tc in clade_buckets[len(rightCladeSet)] if tc.clade ==
                                   rightCladeSet]
                        if matches:
                            for cur_right in matches:
                                # todo this is ugly formatting...
                                # l, r = (left, cur_right) if (
                                #         min(left.clade) < min(cur_right.clade)) \
                                #     else (cur_right, left)
                                l, r = sorted([left, cur_right], key=lambda t: min(t.clade))
                                split = (l, r)
                                if split not in existing_splits:
                                    if compatible(split, clade):
                                        # print(
                                        #     f"-----------------\n"
                                        #     f"Found compatible split:"
                                        #     f"\n{clade}\n"
                                        #     f"{split}\n-----------------"
                                        # )
                                        expanded_clade_partitions[clade].add(split)

    return expanded_clade_partitions


def get_transmission_ccd0(trees, expansion: bool = True):
    m1_observed_clade_counts, m2_observed_clade_split_counts, blockcount_map, branch_lengths_map = (
        get_transmission_maps(trees, type_str="Ancestry"))

    n_trees = len(trees)

    clade_partitions = defaultdict(list)

    for (parent, child1, child2), count in m2_observed_clade_split_counts.items():
        clade_partitions[parent].append((child1, child2))

    if expansion:
        print("Starting expansion")
        expanded_clade_splits = expand_tccd0(clade_partitions)
        clade_partitions = expanded_clade_splits
        print("Expansion has finished...")

    ccd0_probabilities = {}

    def recursive_prob_computer(clade):
        if clade in ccd0_probabilities:
            return ccd0_probabilities[clade]
        clade_value = m1_observed_clade_counts.get(clade, 0) / n_trees

        # leaf
        if len(clade) == 1:
            ccd0_probabilities[clade] = 1.0
            return 1.0

        # cherries are now general clades because leafs can have different transmission ancestors
        transmisison_ccd0_map[clade] = {}
        part_products = []
        for left, right in clade_partitions[clade]:
            # print(f"Recursive work on {clade}")
            left_probability = recursive_prob_computer(left)
            # print(f"Recursive work for left child done...")
            right_probability = recursive_prob_computer(right)
            # print(f"Recursive work for right child done...")
            product = (left_probability * right_probability)
            part_products.append(product)

        total = sum(part_products)

        if total > 0:
            for (left, right), prod in zip(clade_partitions[clade], part_products):
                transmisison_ccd0_map[clade][(left, right)] = prod / total
        else:
            raise ValueError("Need to implement a log fall back...")

        clade_prob = clade_value * total
        ccd0_probabilities[clade] = clade_prob

        child_probabilities = defaultdict(float)

        from collections import deque
        queue = deque()

        # Making sure we only visit each child once, otherwise this explodes, is that correct?...
        # This is different here because now a clade can split differently but with the same clade
        # as one of the children and the other is different, hence we need to check for visited
        # That is. C - C1, C2 but also C - C1, C3 where they have C1 the same but the sibling is
        # different...
        visited = set()

        # Handling multiple possible roots
        for clade in m1_observed_clade_counts.keys():
            if len(clade) == len(trees[0]):
                child_probabilities[clade] = m1_observed_clade_counts[clade] / n_trees
                queue.append(clade)
                visited.add(clade)

        # steps = 0
        while queue:
            # steps += 1
            # if steps % 100 == 0:
            #     print(f"Step {steps}, queue size={len(queue)}")
            #     print(f"Length of stuff: {len(transmisison_ccd0_map.get(clade, {}).values())}")
            clade = queue.popleft()
            parent_prob = child_probabilities[clade]

            for (left, right), ccp in transmisison_ccd0_map.get(clade, {}).items():
                child_probabilities[left] += parent_prob * ccp
                child_probabilities[right] += parent_prob * ccp

                for child in (right, left):
                    if child in transmisison_ccd0_map and child not in visited:
                        # print(f"adding left: {left}")
                        queue.append(child)
                        visited.add(child)

        return clade_prob

    transmisison_ccd0_map = {}
    for clade in m1_observed_clade_counts:
        recursive_prob_computer(clade)
    return transmisison_ccd0_map


def get_map_tree(transmission_ccd0_map):
    seen_resolved_clades = {}
    for cur_clade in sorted(transmission_ccd0_map.keys(), key=len):
        if len(cur_clade) == 2:
            # cherries are now the smallest
            ta_clade_split, prob = max(tccd0_map[cur_clade].items(), key=lambda item: item[1])

            # check if max was unique:
            ties = [k for k, v in tccd0_map[cur_clade].items() if v == prob]
            # if len(ties):
            #     # todo could be annotated here and then if actually used in the MAP tree point it
            #     #  out....
            #     # todo might be valuable to keep all the most likely ones for uncertainty ...
            #     print("TIEBREAKING: Highest probable split of a cherry was not unique")
            seen_resolved_clades[cur_clade] = CladeSplitInfo(ta_clade_split, prob)
        else:
            for cur_clade_split, cur_clade_split_prob in transmission_ccd0_map[cur_clade].items():
                cur_clade_child1, cur_clade_child2 = cur_clade_split
                if len(cur_clade_child1) == 1:
                    cur_clade_child1_prob = 1.0
                else:
                    cur_clade_child1_prob = seen_resolved_clades[cur_clade_child1].prob

                if len(cur_clade_child2) == 1:
                    cur_clade_child2_prob = 1.0
                else:
                    cur_clade_child2_prob = seen_resolved_clades[cur_clade_child2].prob

                cur_clade_split_map_probability = (cur_clade_child1_prob * cur_clade_child2_prob *
                                                   cur_clade_split_prob)
                if cur_clade in seen_resolved_clades:
                    # todo here there can be ties too, maybe need to point those out too...
                    if seen_resolved_clades[cur_clade].prob <= cur_clade_split_map_probability:
                        seen_resolved_clades[cur_clade] = CladeSplitInfo(
                            cur_clade_split, cur_clade_split_map_probability
                        )
                else:
                    seen_resolved_clades[cur_clade] = CladeSplitInfo(
                        cur_clade_split, cur_clade_split_map_probability
                    )

    output = {}
    all_root_clades = [k for k in seen_resolved_clades.keys()
                       if len(k) == len(max(seen_resolved_clades.keys()))]
    max_prob = max(seen_resolved_clades[root].prob for root in all_root_clades)
    best_roots = [root for root in all_root_clades if seen_resolved_clades[root].prob == max_prob]

    if len(best_roots) > 1:
        raise NotImplementedError("There is more than one best root, this needs to be implemented!")

    working_list = [best_roots[0]]
    while working_list:
        cur_parent = working_list.pop()
        split, prob = seen_resolved_clades[cur_parent]
        output[cur_parent] = split
        for split_child in split:
            if len(split_child) > 1:
                working_list.append(split_child)
            else:
                output[split_child] = split_child
    return output


def get_tree_from_dict_of_splits(splits):
    root_split = max(splits.keys())
    output_tree = Tree(support=0, dist=0, name="root")
    output_tree.add_feature("transm_ancest", root_split.transm_ancest)
    icount = 1

    def recursive_children(node, new_split):
        nonlocal splits, icount
        clade_c1, clade_c2 = new_split

        def add_clade(node, clade):
            nonlocal icount, splits
            if len(clade) == 1:
                label = next(iter(clade.clade))
                leaf = node.add_child(name=str(label), dist=1)
                leaf.add_feature("transm_ancest", clade.transm_ancest)
            else:
                internal_node = node.add_child(name=f"child_internal_{icount}", dist=1)
                internal_node.add_feature("transm_ancest", clade.transm_ancest)
                icount += 1
                recursive_children(internal_node, splits[clade])

        add_clade(node, clade_c1)
        add_clade(node, clade_c2)

    recursive_children(output_tree, splits[root_split])

    return output_tree


if __name__ == '__main__':
    from pathlib import Path
    from pyccd import read_nexus_trees

    trees = read_nexus_trees(
        # f"{Path(__file__).parent.absolute().parent.parent}/examples/data/breath32sim.trees",
        # f"{Path(__file__).parent.absolute().parent.parent}/examples/data/roetzer40.trees",
        f"{Path(__file__).parent.absolute().parent.parent}/examples/data/breath32simShort.trees",
        breath_trees=True,
        label_transm_history=True
    )

    # trees = trees[int(len(trees) * 0.95):]
    print(f"Parsed {len(trees)} tree after burnin...")

    tccd0_map = get_transmission_ccd0(trees)

    map_tree = get_map_tree(tccd0_map)

    tree_map_thingy = get_tree_from_dict_of_splits(map_tree)

    print(tree_map_thingy.write(format=5, features=["transm_ancest"], format_root_node=True))
