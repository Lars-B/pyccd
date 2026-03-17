import logging
from collections import defaultdict
from collections import namedtuple

from brokilon.ccd.clades.sranges import SRangesClade

# todo maybe put this with the other in the file: types.py
AncestralSplit = namedtuple(
    "AncestalSplit", ["ancestor", "descendant"]
)


def get_ranges_clade(node, taxon_map):
    ranges_clade = set({})
    for leaf in node:
        if leaf.dist == 0.0:
            if taxon_map[int(leaf.name)].endswith("_first"):
                # assumes that ranges always have the _first (i.e. SAs also do this)
                ranges_clade.add(int(leaf.name))
        else:
            if not taxon_map[int(leaf.name)].endswith("_last"):
                # we ignore all "_last" occurences
                # this is because we want to encode the range as one species, not x_first and x_last
                ranges_clade.add(int(leaf.name))
    return ranges_clade


def get_sranges_map(trees, taxon_map, ccd_type=1):
    clade_count_map = defaultdict(int)
    clade_split_count_map = defaultdict(lambda: defaultdict(int))

    reverse_taxon_map = {value: key for key, value in taxon_map.items()}

    # for node in (node for t in trees for node in t.traverse("levelorder")):
    for t in trees:
        cur_t_last = set()
        for node in t.traverse("levelorder"):

            if len(node.children) == 1:
                logging.debug(f"Found a node with only 1 child")
                continue

            if len(node) > 1:

                # Encoding ranges, int if starting a range with child 1 or 0
                # It will be "End" if we are ending a range with this node
                # None for all the inbetween, i.e. inside and outside ranges

                current_range = None

                if node.children[0].dist == 0.0:
                    if hasattr(node.children[1], "range"):
                        current_range = 1
                    else:
                        current_range = "End"

                if node.children[1].dist == 0.0:
                    if current_range is not None:
                        raise ValueError("This should not happen for my assumption")
                    if hasattr(node.children[0], "range"):
                        current_range = 0
                    else:
                        current_range = "End"

                # only leafs with non zero branch lengths
                if current_range is None:
                    if hasattr(node, "range"):
                        # Using the _first as the range, also using the integer value from nexus file
                        within_range = reverse_taxon_map[f"{node.range}_first"]
                    else:
                        within_range = None

                    c0_leaves = get_ranges_clade(node.children[0], taxon_map)
                    c1_leaves = get_ranges_clade(node.children[1], taxon_map)

                    parent_clade = c0_leaves.union(c1_leaves)
                    if within_range is not None:
                        parent_clade.add(within_range)

                    parent_clade = SRangesClade(
                        frozenset(parent_clade),
                        within_range
                    )

                    clade_count_map[parent_clade] += 1

                    if node.children[0].orientation == "ancestor":
                        assert node.children[1].orientation == "descendant", "Failure...1"

                        anc_clade = SRangesClade(frozenset(c0_leaves), within_range)
                        if within_range:
                            anc_clade = SRangesClade(
                                frozenset(c0_leaves.union({within_range})),
                                within_range
                            )
                        current_split = AncestralSplit(
                            ancestor=anc_clade,
                            descendant=SRangesClade(frozenset(c1_leaves), within_range)
                        )
                    else:
                        assert node.children[1].orientation == "ancestor", "Failure...2"

                        anc_clade = SRangesClade(frozenset(c1_leaves), within_range)
                        if within_range:
                            anc_clade = SRangesClade(
                                frozenset(c1_leaves.union({within_range})),
                                within_range
                            )
                        current_split = AncestralSplit(
                            ancestor=anc_clade,
                            descendant=SRangesClade(frozenset(c0_leaves), within_range)
                        )

                    clade_split_count_map[parent_clade][current_split] += 1

                elif current_range == "End":

                    index_zero_branch = 1 if node.children[1].dist == 0.0 else 0
                    assert node.children[index_zero_branch].dist == 0.0, 'Failure no zero branch.'

                    if not hasattr(node, "range"):
                        # We are at the end of a range, but since the parent
                        # is not inside a range, this is only a sampled ancestor

                        actual_range = int(node.children[index_zero_branch].name)

                        # TODO could use the index zero branch thing here, no need for both
                        non_zero_child_leaves = get_ranges_clade(
                            node.children[0 if index_zero_branch == 1 else 1],
                            taxon_map
                        )

                        parent_clade = SRangesClade(
                            frozenset(non_zero_child_leaves.union({actual_range})),
                            actual_range
                        )
                        anc_clade = SRangesClade(
                            frozenset({}),
                            actual_range
                        )
                        des_clade = SRangesClade(
                            frozenset(non_zero_child_leaves),
                            actual_range
                        )
                        current_split = AncestralSplit(
                            ancestor=anc_clade,
                            descendant=des_clade
                        )

                        clade_count_map[parent_clade] += 1
                        clade_split_count_map[parent_clade][current_split] += 1
                    else:
                        # We end a range and need to deal with adding a proper split
                        # the clade of the parent
                        parent_clade_no_range = get_ranges_clade(node, taxon_map)

                        actual_range = reverse_taxon_map[f"{node.range}_first"]

                        parent_clade = SRangesClade(
                            frozenset(
                                parent_clade_no_range.union(
                                    {actual_range}
                                )
                            ),
                            actual_range
                        )

                        clade_count_map[parent_clade] += 1
                        current_split = AncestralSplit(
                            ancestor=SRangesClade(
                                frozenset({}),
                                actual_range
                            ),
                            descendant=SRangesClade(
                                frozenset(parent_clade_no_range),
                                None
                            )
                        )
                        clade_split_count_map[parent_clade][current_split] += 1

                elif current_range == 0 or current_range == 1:
                    # we are at the start of a range, only have a parent clade at this point...

                    leaves = get_ranges_clade(node.children[current_range], taxon_map)
                    current_ancestral_range = None
                    if hasattr(node.up, "range"):
                        current_ancestral_range = reverse_taxon_map[f"{node.up.range}_first"]

                    leaves.add(reverse_taxon_map[f"{node.children[current_range].range}_first"])

                    current_range_start_clade = SRangesClade(
                        frozenset(leaves),
                        ancestral_range=current_ancestral_range
                    )
                    clade_count_map[current_range_start_clade] += 1
                else:
                    raise ValueError(f"Unexpected value for current range: {current_range}!")

            elif node.is_leaf():

                cur_clade_above = get_ranges_clade(node.up, taxon_map)
                if len(cur_clade_above) == 1:
                    if int(node.name) in cur_clade_above:
                        cur_range = None
                        # Looking above the range if there was an ancestral range
                        if hasattr(node.up.up, "range"):
                            cur_range = node.up.up.range
                            cur_range = reverse_taxon_map[f"{cur_range}_first"]
                        elif node.up.up.children[1].dist == 0:
                            # Sampled ancestor case:
                            cur_range = taxon_map[int(node.up.up.children[1].name)]
                            cur_range = reverse_taxon_map[cur_range.replace("_last", "_first")]
                        elif node.up.up.children[0].dist == 0:
                            cur_range = taxon_map[int(node.up.up.children[0].name)]
                            cur_range = reverse_taxon_map[cur_range.replace("_last", "_first")]
                        leaf_clade = SRangesClade(
                            frozenset({int(node.name)}),
                            ancestral_range=cur_range
                        )
                        clade_count_map[leaf_clade] += 1
                    else:
                        assert taxon_map[int(node.name)].endswith("_last"), "Shouldn't fail?"
                        pass
                else:
                    if taxon_map[int(node.name)].endswith("_last"):
                        assert hasattr(node, "range"), "Needs this!"
                        leaf_clade = SRangesClade(
                            frozenset({reverse_taxon_map[f"{node.range}_first"]}),
                            ancestral_range=reverse_taxon_map[f"{node.range}_first"]

                        )
                        clade_count_map[leaf_clade] += 1
                    else:
                        cur_range = None
                        if hasattr(node.up, "range"):
                            cur_range = reverse_taxon_map[f"{node.up.range}_first"]
                        leaf_clade = SRangesClade(
                            frozenset({int(node.name)}),
                            ancestral_range=cur_range
                        )
                        clade_count_map[leaf_clade] += 1

    # Return a regular dict to avoid adding things later that shouldn't be in there.
    return dict(clade_count_map), {k: dict(v) for k, v in clade_split_count_map.items()}


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

        # todo we now get current_clade that have no splits associated with them,
        #  fix the clade extraction function further, leaves should be fixed now.
        #  One of the problem cases: taxon_map[81], taxon_map[83], both _first of differnt species,
        #  No ancestral range, might be causing it because there should be a range? another SA case?
        for current_split in clade_split_count_map[current_clade]:

            anc_prob, desc_prob = 0, 0

            if len(current_split.ancestor.clade) == 1:
                leaf_observations = sum(
                    [clade_count_map[k] for k in leaf_clades if
                     k.clade == current_split.ancestor.clade]
                )
                anc_prob = clade_count_map[current_split.ancestor] / leaf_observations
                assert anc_prob <= 1.0, "Prob failure1..."
            elif len(current_split.ancestor.clade) == 0:
                anc_prob = 1
            else:
                anc_prob = seen_resolved_clades[current_split.ancestor][0]

            if len(current_split.descendant.clade) == 1:
                leaf_observations = sum(
                    [clade_count_map[k] for k in leaf_clades if
                     k.clade == current_split.descendant.clade]
                )
                desc_prob = clade_count_map[current_split.descendant] / leaf_observations
                assert desc_prob <= 1.0, "Prob failure2..."
            elif len(current_split.descendant.clade) == 0:
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
                        chose_prob, chosen_split = seen_resolved_clades[current_clade][:2]

                    # resolving tiebreak in the seen_resolved_clades, True indicates the tiebreak
                    seen_resolved_clades[current_clade] = (chose_prob, chosen_split, True)
            else:
                seen_resolved_clades[current_clade] = (split_probability, current_split, False)
    # End of seen_resolved_clades construction

    print("help me.")

    return None

# todo we will need to do a second post processing step over all ranges
#  that ended up being leaves without first and last, there add ranges only at the leaves
# todo maybe keep a dict of int --> is_range? or something like that?
