from collections import defaultdict

from brokilon.ccd.clades.sranges import SRangesClade
from brokilon.core import Tree

import inspect
from dataclasses import dataclass, field


def get_callsite():
    frame = inspect.stack()[2]  # caller of caller
    return f"{frame.filename}:{frame.lineno} ({frame.function})"


@dataclass(frozen=True)
class AncestralSplit:
    ancestor: SRangesClade
    descendant: SRangesClade
    source: str = field(default_factory=get_callsite)


def prelabel_tree(tree, taxon_map):
    opened_ranges = set({})
    closed_ranges = set({})
    sampled_ancestors = set({})

    # This requires a two pass thing because the range can be opend and closed
    # at the same distance (in terms of number of nodes/edges) to the root
    # This could be avoided by doing a levelorder algorithm
    # that also takes branch lengths into account
    # Or first sorting the leafs by distance to root and then iterating over that list instead.
    for leaf in tree:
        current_taxon = taxon_map[int(leaf.name)]
        if leaf.dist == 0.0:
            # non regular leaf
            if current_taxon.endswith("_first"):
                other_leaves = [taxon_map[int(t.name)] for t in leaf.up if t is not leaf]

                range_end = current_taxon.replace("_first", "_last")

                if range_end in other_leaves:
                    # We are at the start of a range
                    opened_ranges.add(current_taxon)

                    if not hasattr(leaf.up, "rangetype"):
                        leaf.up.add_feature("rangetype", "range_start")
                    else:
                        assert leaf.up.rangetype == "range_end", "Failure"
                        leaf.up.rangetype = "leaf_range"
                    leaf.up.add_feature("range", current_taxon[:-6])
                    leaf.add_feature("rangetype", "range_start")
                else:
                    # We are looking at a SA
                    leaf.add_feature("rangetype", "sampled_ancestor")
                    leaf.up.add_feature("rangetype", "sampled_ancestor")

                    assert not hasattr(leaf.up, "range"), "This should be SA!"
                    # Adding the SA as an attribute to the parent for further use
                    leaf.up.add_feature("sampledancestor", taxon_map[int(leaf.name)][:-6])

                    sampled_ancestors.add(current_taxon)
            elif current_taxon.endswith("_last"):
                # This should always be the end of a range

                if not hasattr(leaf.up, "rangetype"):
                    leaf.up.add_feature("rangetype", "range_end")
                else:
                    assert leaf.up.rangetype == "range_start", "failure"
                    leaf.up.rangetype = "leaf_range"
                leaf.add_feature("rangetype", "range_end")

                closed_ranges.add(current_taxon)
            else:
                raise ValueError("This should not happen")
    # second pass to catch all the closed ranges that happen at the same node as opening...
    for leaf in tree:
        current_taxon = taxon_map[int(leaf.name)]
        if current_taxon.endswith("_last"):
            # Closing ranges with leafs that don't have dist 0
            if current_taxon.replace("_last", "_first") in opened_ranges:
                if not hasattr(leaf.up, "rangetype"):
                    if leaf.dist == 0.0:
                        leaf.up.add_feature("rangetype", "range_end")
                elif leaf.up.rangetype == "range_start":
                    # The only case that we need to chnage it to leaf_range is this:
                    # It is already a range_start, and we would overwrite it with range_end
                    leaf.up.rangetype = "leaf_range"

                leaf.add_feature("rangetype", "range_end")

                closed_ranges.add(current_taxon)

            else:
                # assuming that leaf is not end of a range.
                continue
    if len(opened_ranges) != len(closed_ranges):
        raise ValueError("This is not possible and a bug is found.")
    print(f"Tree has {len(sampled_ancestors)} SAs and {len(opened_ranges)} ranges.")


def get_sranges_map(trees, taxon_map, ccd_type=1):
    clade_count_map = defaultdict(int)
    clade_split_count_map = defaultdict(lambda: defaultdict(int))

    reverse_taxon_map = {value: key for key, value in taxon_map.items()}

    # for node in (node for t in trees for node in t.traverse("levelorder")):
    for t in trees:
        prelabel_tree(t, taxon_map)

        for node in t.traverse("levelorder"):
            # We can ignore all leafs that are 0.0, i.e. internal leafs that encode a non leaf...
            if node.is_leaf():
                if node.dist == 0.0:
                    continue
                if hasattr(node, "rangetype"):
                    if node.rangetype == "range_end":
                        if hasattr(node.up, "rangetype"):
                            if node.up.rangetype == "leaf_range":
                                continue
                        assert node.dist != 0.0, "If this happens we need to not allow it?!"
                        assert node.up.range == node.range, "The ranges don't match?"
                        assert (taxon_map[int(node.name)].replace("_last", "_first") ==
                                f"{node.range}_first"), \
                            "More problems"
                        clade_count_map[SRangesClade(
                            frozenset({f"{node.range}_first"}),
                            f"{node.range}_first"
                        )] += 1
                    else:
                        # ignore range_start and sampled ancestors
                        continue
                else:
                    # regular leaf
                    assert taxon_map[int(node.name)].endswith("_first"), "Failed leaf case..."
                    cur_range = None
                    # todo double check this with the following cases, is there some redundancy
                    #  comparing to the match case below?
                    if hasattr(node.up, "range"):
                        cur_range = f"{node.up.range}_first"
                    elif hasattr(node.up, "rangetype"):
                        match node.up.rangetype:
                            case "range_end":
                                cur_range = f"{node.up.range}_first"
                            case "sampled_ancestor":
                                cur_range = f"{node.up.sampledancestor}_first"
                            case _:
                                raise AssertionError(
                                    "If this happens there might be a missing case..."
                                )

                    new_clade = SRangesClade(
                        frozenset({taxon_map[int(node.name)]}),
                        cur_range
                    )
                    clade_count_map[new_clade] += 1

            elif hasattr(node, "rangetype"):
                match node.rangetype:
                    case "sampled_ancestor":
                        current_range = None
                        if hasattr(node.up, "range"):
                            current_range = f"{node.up.range}_first"
                        elif hasattr(node.up, "sampledancestor"):
                            current_range = f"{node.up.sampledancestor}_first"
                        parent_clade_set = {taxon_map[int(l.name)].replace("_last", "_first") for l
                                            in node}

                        child0_clade_set = {taxon_map[int(l.name)].replace("_last", "_first") for l
                                            in node.children[0]}
                        child1_clade_set = {taxon_map[int(l.name)].replace("_last", "_first") for l
                                            in node.children[1]}

                        assert parent_clade_set == child0_clade_set.union(child1_clade_set)

                        current_parent_clade = SRangesClade(
                            frozenset(parent_clade_set),
                            current_range
                        )

                        if node.children[0].dist == 0.0:
                            current_split = AncestralSplit(
                                ancestor=SRangesClade(
                                    frozenset({}),
                                    taxon_map[int(node.children[0].name)]
                                ),
                                descendant=SRangesClade(
                                    frozenset(child1_clade_set),
                                    taxon_map[int(node.children[0].name)]
                                )
                            )
                        elif node.children[1].dist == 0.0:
                            current_split = AncestralSplit(
                                ancestor=SRangesClade(
                                    frozenset({}),
                                    taxon_map[int(node.children[1].name)]
                                ),
                                descendant=SRangesClade(
                                    frozenset(child0_clade_set),
                                    taxon_map[int(node.children[1].name)]
                                )
                            )

                        clade_count_map[current_parent_clade] += 1
                        clade_split_count_map[current_parent_clade][current_split] += 1
                    case "range_start":
                        current_range = None
                        # if hasattr(node.up, "range"):
                        #     current_range = f"{node.up.range}_first"
                        # elif hasattr(node.up, "sampledancestor"):
                        #     current_range = f"{node.up.sampledancestor}_first"

                        if node.up:
                            if hasattr(node.up, "rangetype"):
                                match node.up.rangetype:
                                    case "range_end":
                                        current_range = f"{node.up.range}_first"
                                    case "sampled_ancestor":
                                        assert hasattr(node.up, "sampledancestor"),\
                                            "This needs fixing"
                                        current_range = f"{node.up.sampledancestor}_first"
                                    case "range_start":
                                        assert node.up.range == node.range, \
                                            "something wrong with node.up.up case"
                                        # todo check this case again...
                                        if hasattr(node.up.up, "range"):
                                            current_range = f"{node.up.up.range}_first"
                                    case _:
                                        pass
                            elif hasattr(node.up, "range"):
                                current_range = f"{node.up.range}_first"

                        parent_clade = SRangesClade(
                            frozenset(
                                {taxon_map[int(l.name)].replace("_last", "_first") for l in node}
                            ),
                            current_range
                        )
                        clade_count_map[parent_clade] += 1
                    case "range_end":
                        # Keeping just the taxon name as range for now...
                        child0_clade_set = {taxon_map[int(l.name)].replace("_last", "_first") for l
                                            in node.children[0]}
                        child1_clade_set = {taxon_map[int(l.name)].replace("_last", "_first") for l
                                            in node.children[1]}
                        parent_clade_set = {
                            taxon_map[int(l.name)].replace("_last", "_first") for l in node
                        }

                        assert parent_clade_set == child0_clade_set.union(child1_clade_set), \
                            "Failure in range_end case..."

                        # The following is a range end, hence one of the two clade sets
                        # is the one taxa that is coding the end of the range ..._last
                        # We need to figure out which one, set it as descendant,
                        # and then we remove it from the taxon set to encode the end of range

                        # if node.children[0].orientation == "ancestor":
                        if child1_clade_set == {f"{node.range}_first"}:
                            current_split = AncestralSplit(
                                ancestor=SRangesClade(
                                    frozenset(child0_clade_set),
                                    f"{node.range}_first"
                                ),
                                descendant=SRangesClade(
                                    # frozenset(child1_clade_set),
                                    frozenset({}),
                                    f"{node.range}_first"
                                )
                            )
                        elif child0_clade_set == {f"{node.range}_first"}:
                            current_split = AncestralSplit(
                                ancestor=SRangesClade(
                                    frozenset(child1_clade_set),
                                    f"{node.range}_first"
                                ),
                                descendant=SRangesClade(
                                    # frozenset(child0_clade_set),
                                    frozenset({}),
                                    f"{node.range}_first"
                                )
                            )
                        else:
                            raise AssertionError(f"Failure: range end should have "
                                                 f"single taxon leaf set which is the range...")
                        parent_range = None
                        if node.up:
                            if hasattr(node.up, "rangetype"):
                                match node.up.rangetype:
                                    # todo this case with the the case elif below? redundancy?
                                    case "range_end":
                                        parent_range = f"{node.up.range}_first"
                                    case "sampled_ancestor":
                                        assert hasattr(node.up, "sampledancestor"), \
                                            "This needs fixing"
                                        parent_range = f"{node.up.sampledancestor}_first"
                                    case "range_start":
                                        assert node.up.range == node.range, \
                                            "something wrong with node.up.up case"
                                        if hasattr(node.up.up, "range"):
                                            parent_range = f"{node.up.up.range}_first"
                                        elif hasattr(node.up.up, "sampledancestor"):
                                            parent_range = f"{node.up.up.sampledancestor}_first"
                                    case _:
                                        pass
                            elif hasattr(node.up, "range"):
                                parent_range = f"{node.up.range}_first"

                        parent_clade = SRangesClade(
                            frozenset(parent_clade_set),
                            parent_range
                        )
                        clade_count_map[parent_clade] += 1
                        clade_split_count_map[parent_clade][current_split] += 1
                    case "leaf_range":
                        # We are looking at a node that is parent of a range that is a leaf
                        # We don't need to add a split here?
                        # We only count the clade, which is the taxon and itself as a range?

                        # todo fix this to something less messy...
                        current_leaf = [l.name for l in node if l.dist == 0.0][0]
                        range_taxon = reverse_taxon_map[f"{node.range}_first"]

                        # todo what is this assert doing?
                        assert range_taxon == int(current_leaf), \
                            (f"Failure: {range_taxon} != {current_leaf}")

                        # Leaf ranges can not have themselves as a range,
                        # but there can be other ranges as ancestors above
                        # Leaf ranges have no ancestral range of themselves because
                        # otherwise this information is assumed to be present above the leaf
                        # and therefore breaks the logic
                        # One way to avoid this and fully encode leaf ranges is
                        # to have _first and _last taxon as part of the clades...
                        current_range = None
                        if hasattr(node.up, "range"):
                            current_range = f"{node.up.range}_first"
                        elif hasattr(node.up, "sampledancestor"):
                            current_range = f"{node.up.sampledancestor}_first"
                        clade_count_map[
                            SRangesClade(
                                frozenset({taxon_map[int(current_leaf)]}),
                                current_range
                            )
                        ] += 1
                    case _:
                        raise ValueError("Unsupported value...")

            else:
                parent_clade_set = {
                    taxon_map[int(l.name)].replace("_last", "_first") for l in node
                }

                child0_clade_set = {
                    taxon_map[int(l.name)].replace("_last", "_first") for l in node.children[0]
                }
                child1_clade_set = {
                    taxon_map[int(l.name)].replace("_last", "_first") for l in node.children[1]
                }

                assert child1_clade_set.union(child0_clade_set) == parent_clade_set, \
                    "Something is wrong with the clades..."
                cur_range = None
                parent_range = None
                if hasattr(node, "range"):
                    cur_range = f"{node.range}_first"
                if node.up:
                    if hasattr(node.up, "rangetype"):
                        match node.up.rangetype:
                            case "range_end":
                                parent_range = f"{node.up.range}_first"
                            case "sampled_ancestor":
                                assert hasattr(node.up, "sampledancestor"), "This needs fixing"
                                parent_range = f"{node.up.sampledancestor}_first"
                            case "range_start":
                                assert node.up.range == node.range, "something wrong with node.up.up case"
                                if hasattr(node.up.up, "range"):
                                    parent_range = f"{node.up.up.range}_first"
                                elif hasattr(node.up.up, "sampledancestor"):
                                    parent_range = f"{node.up.up.sampledancestor}_first"
                            case _:
                                pass
                    elif hasattr(node.up, "range"):
                        parent_range = f"{node.up.range}_first"

                parent_clade = SRangesClade(frozenset(parent_clade_set), parent_range)

                if node.children[0].orientation == "ancestor":
                    current_split = AncestralSplit(
                        ancestor=SRangesClade(frozenset(child0_clade_set), cur_range),
                        descendant=SRangesClade(frozenset(child1_clade_set), cur_range),
                    )
                else:
                    current_split = AncestralSplit(
                        ancestor=SRangesClade(frozenset(child1_clade_set), cur_range),
                        descendant=SRangesClade(frozenset(child0_clade_set), cur_range),
                    )

                clade_count_map[parent_clade] += 1
                clade_split_count_map[parent_clade][current_split] += 1

    return (dict(clade_count_map),
            {clade: dict(splits) for clade, splits in clade_split_count_map.items()})


def get_sranges_map_tree(
        clade_count_map,
        clade_split_count_map,
        taxon_map, reverse_taxon_map
):
    seen_resolved_clades = {}
    leaf_clades = []
    for current_clade in sorted(clade_count_map.keys(), key=len):
        if len(current_clade) == 1:
            leaf_clades.append(current_clade)
            continue
        elif len(current_clade) == 0:
            # this is splitting off a range or SA
            continue

        for current_split in clade_split_count_map[current_clade]:

            anc_prob, desc_prob = 0, 0

            if len(current_split.ancestor.clade) == 1:
                leaf_observations = sum(
                    [clade_count_map[k] for k in leaf_clades if
                     k.clade == current_split.ancestor.clade]
                )

                anc_prob = clade_count_map[current_split.ancestor] / leaf_observations
                assert 0 <= anc_prob <= 1.0, "Prob failure1..."
            elif len(current_split.ancestor.clade) == 0:
                anc_prob = 1
            else:
                # [c for c in seen_resolved_clades if c.clade == current_split.ancestor.clade], current_split.ancestor, seen_resolved_clades[current_split.ancestor]
                anc_prob = seen_resolved_clades[current_split.ancestor][0]

            if len(current_split.descendant.clade) == 1:
                leaf_observations = sum(
                    [clade_count_map[k] for k in leaf_clades if
                     k.clade == current_split.descendant.clade]
                )

                # current_split.descendant, [k for k in clade_count_map if k.clade == current_split.descendant.clade]
                desc_prob = clade_count_map[current_split.descendant] / leaf_observations
                assert 0 <= desc_prob <= 1.0, "Prob failure2..."
            elif not current_split.descendant.clade:
                desc_prob = 1
            else:
                desc_prob = seen_resolved_clades[current_split.descendant][0]

            split_probability = anc_prob * desc_prob * (
                    clade_split_count_map[current_clade][current_split] /
                    clade_count_map[current_clade]
            )

            if current_clade in seen_resolved_clades:
                if seen_resolved_clades[current_clade][0] < split_probability:
                    seen_resolved_clades[current_clade] = (split_probability, current_split, False)
                elif seen_resolved_clades[current_clade][0] == split_probability:
                    # Tie breaking randomly, either keep the old or pick the new
                    import random
                    if random.random() < 0.5:
                        chosen_prob, chosen_split = split_probability, current_split
                    else:
                        chosen_prob, chosen_split = seen_resolved_clades[current_clade][:2]

                    # resolving tiebreak in the seen_resolved_clades, True indicates the tiebreak
                    seen_resolved_clades[current_clade] = (chosen_prob, chosen_split, True)
            else:
                seen_resolved_clades[current_clade] = (split_probability, current_split, False)
    # End of seen_resolved_clades construction

    # todo we should use the named tuple to store prob and splits etc...

    output = {}
    max_key_value = len(max(seen_resolved_clades.keys()))
    all_root_clades = [k for k in seen_resolved_clades if len(k) == max_key_value]
    max_root_prob = max(seen_resolved_clades[root][0] for root in all_root_clades)
    best_roots = [root for root in all_root_clades
                  if seen_resolved_clades[root][0] == max_root_prob]

    if len(best_roots) > 1:
        raise NotImplementedError("Currently no tie breaking support for multiple roots")
    else:
        root = best_roots[0]

    map_tree = get_sranges_tree_from_seen_resolved_clades(seen_resolved_clades, root)

    return map_tree


def get_sranges_tree_from_seen_resolved_clades(
        seen_resolved_clades,
        root_clade
):
    # todo ignoring tie breaking for now...
    # todo ignoring ranges at the moment...

    out_tree = Tree(
        support=100,
        dist=1.0,
        name="root"
    )

    icount = 1

    def recursive_sranges_children(node, split):
        nonlocal icount, seen_resolved_clades

        cur_anc = split.ancestor
        cur_descendant = split.descendant

        def add_clade(parent_node, clade, orientation):
            nonlocal icount, seen_resolved_clades

            # leaf
            if len(clade.clade) == 1:
                label = next(iter(clade.clade))

                leaf = parent_node.add_child(
                    name=label,
                    dist=10,
                    # mean(branch_length_map.get(clade, [1])),
                    support=1.0,
                )
                # add feature ancestor or descendant
                leaf.add_feature("orientation", orientation)
                leaf.add_feature('ancestral_range', clade.ancestral_range)
            # internal node
            else:
                internal_node = parent_node.add_child(
                    name=f"internal_{icount}",
                    dist=10,
                    # mean(branch_length_map.get(clade, [1])),
                    support=10
                    # split_support_map.get(split, 1.0),
                )

                internal_node.add_feature("orientation", orientation)
                internal_node.add_feature("ancestral_range", clade.ancestral_range)

                icount += 1

                # recurse if this clade was resolved
                if clade in seen_resolved_clades:
                    child_split = seen_resolved_clades[clade][1]
                    recursive_sranges_children(internal_node, child_split)

        add_clade(node, cur_anc, "ancestor")
        add_clade(node, cur_descendant, "descendant")

    recursive_sranges_children(
        out_tree,
        seen_resolved_clades[root_clade][1]
    )

    return out_tree
