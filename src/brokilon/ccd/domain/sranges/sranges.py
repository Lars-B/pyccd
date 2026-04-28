from collections import defaultdict
from collections import namedtuple

from brokilon.ccd.clades.sranges import SRangesClade

# todo maybe put this with the other in the file: types.py
AncestralSplit = namedtuple(
    "AncestalSplit", ["ancestor", "descendant"]
)


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
                    assert leaf.up.rangetype == "range_start", "failue"
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
            if f"{current_taxon[:-5]}_first" in opened_ranges:
                if not hasattr(leaf.up, "rangetype"):
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

                            # todo add more cases, i think the problem comes from SAs so that
                            #  is the case I need to fix...

                        assert node.dist != 0.0, "If this happens we need to not allow it?!"
                        assert node.up.range == node.range, "The ranges don't match?"
                        assert (taxon_map[int(node.name)].replace("_last", "_first") ==
                                f"{node.range}_first"), \
                            "More problems"
                        clade_count_map[SRangesClade(
                            frozenset({f"{node.range}_first"}),
                            node.range
                        )] += 1
                    else:
                        # ignore range_start and sampled ancestors
                        continue
                else:
                    # regular leaf
                    assert taxon_map[int(node.name)].endswith("_first"), "Failed leaf case..."
                    clade_count_map[SRangesClade(
                        frozenset({taxon_map[int(node.name)]}),
                        None
                    )] += 1

            elif hasattr(node, "rangetype"):
                match node.rangetype:
                    case "sampled_ancestor":
                        current_range = None
                        if hasattr(node.up, "range"):
                            current_range = node.up.range
                        parent_clade_set = {taxon_map[int(l.name)].replace("_last", "_first") for l in node}

                        child0_clade_set = {taxon_map[int(l.name)].replace("_last", "_first") for l in node.children[0]}
                        child1_clade_set = {taxon_map[int(l.name)].replace("_last", "_first") for l in node.children[1]}

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
                        parent_clade = SRangesClade(
                            frozenset(
                                {taxon_map[int(l.name)].replace("_last", "_first") for l in node}
                            ),
                            node.range
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

                        if node.children[0].orientation == "ancestor":
                            current_split = AncestralSplit(
                                    ancestor=SRangesClade(
                                        frozenset(child0_clade_set),
                                        node.range
                                    ),
                                    descendant=SRangesClade(
                                        frozenset(child1_clade_set),
                                        node.range
                                    )
                                )
                        else:
                            current_split = AncestralSplit(
                                    ancestor=SRangesClade(
                                        frozenset(child1_clade_set),
                                        node.range
                                    ),
                                    descendant=SRangesClade(
                                        frozenset(child0_clade_set),
                                        node.range
                                    )
                                )
                        clade_split_count_map[
                            SRangesClade(
                                frozenset(parent_clade_set),
                                node.range)
                        ][current_split] += 1
                    case "leaf_range":
                        # We are looking at a node that is parent of a range that is a leaf
                        # We don't need to add a split here?
                        # We only count the clade, which is the taxon and itself as a range?

                        # This can surely be done in a nicer way...
                        current_leaf = [l.name for l in node if l.dist == 0.0][0]
                        range_taxon = reverse_taxon_map[f"{node.range}_first"]
                        assert range_taxon == int(current_leaf), (f"Failure: "
                                                                  f"{range_taxon} != {current_leaf}")

                        # Leaf ranges have no ancestral range, should be encoded in leaf is a range
                        # otherwise this information is not present above the leaf
                        # and therefore breaks the logic
                        clade_count_map[
                            SRangesClade(
                                frozenset({taxon_map[int(current_leaf)]}),
                                # node.range
                                None
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
                if hasattr(node, "range"):
                    cur_range = node.range
                else:
                    if node.up:
                        if hasattr(node.up, "rangetype"):
                            if node.up.rangetype == "range_end":
                                cur_range = node.up.range
                            elif node.up.rangetype == "sampled_ancestor":
                                assert hasattr(node.up, "sampledancestor"), "This needs fixing"
                                cur_range = node.up.sampledancestor
                parent_clade = SRangesClade(frozenset(parent_clade_set), cur_range)

                if node.children[0].orientation == "ancestor":
                    cur_split = AncestralSplit(
                        ancestor=SRangesClade(frozenset(child0_clade_set), cur_range),
                        descendant=SRangesClade(frozenset(child1_clade_set), cur_range),
                    )
                else:
                    cur_split = AncestralSplit(
                        ancestor=SRangesClade(frozenset(child1_clade_set), cur_range),
                        descendant=SRangesClade(frozenset(child0_clade_set), cur_range),
                    )

                clade_count_map[parent_clade] += 1
                clade_split_count_map[parent_clade][cur_split] += 1

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

                # Todo this is still wrong i believe, we are missing that the leaf without anything
                #  could still be counted sometimes, needs to be fixed but this is a quick debug fix
                if leaf_observations == 0:
                    anc_prob = 1
                    raise ValueError("I don't think this is necessary anymore?")
                else:
                    anc_prob = clade_count_map[current_split.ancestor] / leaf_observations
                assert anc_prob <= 1.0, "Prob failure1..."
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
                # todo this is wrong, see above...
                if leaf_observations == 0:
                    desc_prob = 1
                    raise ValueError("I don't think this is necessary anymore?")
                else:
                    # todo clade count map contains clade without range whereas the split has the range
                    #  this seems like an error that is fixable,
                    #  might have some logic problem when adding leaf clades to clade count map?
                    #  Probably add some debug info into the clades, i.e. another attribute like when it was added to the clade_count_map to find where things go wrong...

                    # current_split.descendant, [k for k in clade_count_map if k.clade == current_split.descendant.clade]
                    desc_prob = clade_count_map[current_split.descendant] / leaf_observations
                assert desc_prob <= 1.0, "Prob failure2..."
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

    return None
