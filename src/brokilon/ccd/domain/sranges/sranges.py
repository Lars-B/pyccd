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
            if taxon_map[int(leaf.name)].endswith("_last"):
                # if we have last we add the one with first
                first_taxon = next(k for k, v in taxon_map.items()
                                   if v == taxon_map[int(leaf.name)].replace("_last", "_first"))
                ranges_clade.add(first_taxon)
            else:
                # if we don't have _last just add it
                ranges_clade.add(int(leaf.name))
    return ranges_clade


def get_sranges_map(trees, taxon_map, ccd_type=1):
    clade_count_map = defaultdict(int)
    clade_split_count_map = defaultdict(lambda: defaultdict(int))

    reverse_taxon_map = {value: key for key, value in taxon_map.items()}

    # for node in (node for t in trees for node in t.traverse("levelorder")):
    for t in trees:
        cur_leaf_list = [l.name for l in t]

        # Could make a dict for each tree, for more debugging stuff...
        observed_clades = set()

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
                    current_node_range = None
                    if hasattr(node, "range"):
                        # Using the _first as the range, also using the integer value from nexus file
                        current_node_range = reverse_taxon_map[f"{node.range}_first"]

                    c0_leaves = get_ranges_clade(node.children[0], taxon_map)
                    c1_leaves = get_ranges_clade(node.children[1], taxon_map)

                    parent_clade_set = c0_leaves.union(c1_leaves)
                    ancestral_range = current_node_range
                    if current_node_range:
                        # this is if the current node is inside a range
                        parent_clade_set.add(current_node_range)
                    elif not node.is_root():
                        if hasattr(node.up, "range"):
                            # the node above the current one was within a range i.e. it ends a range
                            ancestral_range = reverse_taxon_map[f"{node.up.range}_first"]
                        elif node.up.children[0].dist == 0.0:
                            # the above node had child 0 which is a sampled ancestor
                            assert node.up.children[0].is_leaf(), "Failure SA1"
                            assert taxon_map[int(node.up.children[0].name)].endswith("_first")
                            ancestral_range = int(node.up.children[0].name)
                        elif node.up.children[1].dist == 0.0:
                            # the above node had child 1 which is a sampled ancestor
                            assert node.up.children[1].is_leaf(), "Failure SA2"
                            assert taxon_map[int(node.up.children[1].name)].endswith("_first")
                            ancestral_range = int(node.up.children[1].name)

                    parent_clade = SRangesClade(
                        frozenset(parent_clade_set),
                        ancestral_range
                    )

                    clade_count_map[parent_clade] += 1
                    observed_clades.add(parent_clade)

                    if node.children[0].orientation == "ancestor":
                        assert node.children[1].orientation == "descendant", "Failure...1"

                        anc_clade = SRangesClade(frozenset(c0_leaves), current_node_range)

                        current_split = AncestralSplit(
                            ancestor=anc_clade,
                            descendant=SRangesClade(frozenset(c1_leaves), current_node_range)
                        )
                    else:
                        assert node.children[1].orientation == "ancestor", "Failure...2"

                        anc_clade = SRangesClade(frozenset(c1_leaves), current_node_range)

                        current_split = AncestralSplit(
                            ancestor=anc_clade,
                            descendant=SRangesClade(frozenset(c0_leaves), current_node_range)
                        )

                    clade_split_count_map[parent_clade][current_split] += 1

                elif current_range == "End":

                    index_zero_branch = 1 if node.children[1].dist == 0.0 else 0
                    assert node.children[index_zero_branch].dist == 0.0, 'Failure no zero branch.'

                    if not hasattr(node, "range"):
                        # We are at the end of a range, but since the parent
                        # is not inside a range, this is only a sampled ancestor

                        sampled_ancestor = int(node.children[index_zero_branch].name)
                        clade_at_node = get_ranges_clade(node, taxon_map)

                        node_ancestral_range = None
                        if hasattr(node.up, "range"):
                            node_ancestral_range = reverse_taxon_map[f"{node.up.range}_first"]

                        parent_clade = SRangesClade(
                            frozenset(clade_at_node),
                            node_ancestral_range
                        )
                        anc_clade = SRangesClade(
                            frozenset({}),
                            sampled_ancestor
                        )
                        des_clade = SRangesClade(
                            frozenset(clade_at_node - {sampled_ancestor}),
                            sampled_ancestor
                        )
                        current_split = AncestralSplit(
                            ancestor=anc_clade,
                            descendant=des_clade
                        )

                        clade_count_map[parent_clade] += 1
                        observed_clades.add(parent_clade)
                        clade_split_count_map[parent_clade][current_split] += 1
                    else:
                        # We end a range and need to deal with adding a proper split
                        # the clade of the parent
                        parent_clade_no_range = get_ranges_clade(node, taxon_map)

                        # Clarification of the following up and up.up
                        # We know that we are at a node that ends a range
                        # that means one child has dist == 0.0
                        # so our split will have one empty clade with this range as ancestor
                        # Furthermore, we need to find ancestor for current node
                        # Either the node above was in a range, which will then be our current range
                        # If the node above was not in a range then
                        #   that node started the range ending at node
                        # If that is the case, our split is for the clade that is below node.up
                        # and the ancestor for the clade at node.up is coming from the range of node.up.up
                        # TODO this is missing SA cases at node.up.up which might cause problems...
                        #  add more if, SA if not range but one node is _first with dist 0.0????
                        # End clarification

                        # default no ancestral range
                        actual_range = None
                        if hasattr(node.up, "range"):
                            # either the node above is part of the range and has some split
                            actual_range = reverse_taxon_map[f"{node.up.range}_first"]
                        elif hasattr(node.up.up, "range"):
                            # or the node above is the start of the range, hence the range should be
                            # encoded by up.up.range
                            # todo this might miss sampled ancestors?
                            # TODO REALLY?
                            actual_range = reverse_taxon_map[f"{node.up.up.range}_first"]
                        ending_range = reverse_taxon_map[f"{node.range}_first"]
                        parent_clade = SRangesClade(
                            frozenset(
                                parent_clade_no_range.union({ending_range})
                            ),
                            actual_range
                        )

                        # todo is this double counting something? I think so...
                        # if parent_clade not in observed_clades:
                        #     clade_count_map[parent_clade] += 1
                        #     observed_clades.add(parent_clade)
                        # else:
                        #     print("Double counting stuff...")

                        # this is the range that currently ends
                        desc = SRangesClade(
                                frozenset(parent_clade_no_range),
                                ending_range
                            )
                        current_split = AncestralSplit(
                            ancestor=SRangesClade(
                                frozenset({}),
                                ending_range
                            ),
                            descendant=desc
                        )
                        clade_count_map[desc] += 1
                        observed_clades.add(desc)
                        clade_split_count_map[parent_clade][current_split] += 1

                elif current_range == 0 or current_range == 1:
                    # we can also have this when we are looking at sampled ancestors
                    range_start_node = 1 if current_range == 0 else 0
                    assert int(node.children[range_start_node].name) in taxon_map
                    range_start_name = taxon_map[int(node.children[range_start_node].name)]
                    assert range_start_name.endswith("_first")

                    leaves = get_ranges_clade(node.children[current_range], taxon_map)

                    if not (reverse_taxon_map[f"{range_start_name[:-5]}last"]
                            in cur_leaf_list) and len(leaves) > 1:

                        current_sampled_ancstor = reverse_taxon_map[range_start_name]
                        # The sampled ancestor is part of the clade at the current clade, but splits off
                        current_ancestral_range = None
                        if hasattr(node.up, "range"):
                            current_ancestral_range = reverse_taxon_map[f"{node.up.range}_first"]
                        current_clade_with_sa = SRangesClade(
                            frozenset(leaves.union({current_sampled_ancstor})),
                            ancestral_range=current_ancestral_range
                        )

                        clade_count_map[current_clade_with_sa] += 1
                        observed_clades.add(current_clade_with_sa)

                        # Since we have only a SA we need to also record the split that happens
                        clade_after_sa_split = SRangesClade(
                            frozenset(leaves - {current_sampled_ancstor}),
                            ancestral_range=current_sampled_ancstor
                        )

                        clade_count_map[clade_after_sa_split] += 1
                        observed_clades.add(clade_after_sa_split)

                        current_split = AncestralSplit(
                            ancestor=clade_after_sa_split,
                            descendant=SRangesClade(None, current_sampled_ancstor)
                        )
                        clade_split_count_map[current_clade_with_sa][current_split] += 1
                    else:
                        # This is the start of an actual range
                        # we are at the start of a range, only have a parent clade at this point...

                        # leaves = get_ranges_clade(node.children[current_range], taxon_map)
                        current_ancestral_range = None
                        if hasattr(node.up, "range"):
                            current_ancestral_range = reverse_taxon_map[f"{node.up.range}_first"]

                        leaves.add(reverse_taxon_map[f"{node.children[current_range].range}_first"])

                        current_range_start_clade = SRangesClade(
                            frozenset(leaves),
                            ancestral_range=current_ancestral_range
                        )
                        clade_count_map[current_range_start_clade] += 1
                        observed_clades.add(current_range_start_clade)
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
                        observed_clades.add(leaf_clade)
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
                        observed_clades.add(leaf_clade)
                    else:
                        cur_range = None
                        if hasattr(node.up, "range"):
                            cur_range = reverse_taxon_map[f"{node.up.range}_first"]
                        elif node.up.children[0].dist == 0:
                            assert taxon_map[int(node.up.children[0].name)].endswith("_first"), \
                                "SA1 failure"
                            cur_range = int(node.up.children[0].name)

                        elif node.up.children[1].dist == 0:
                            assert taxon_map[int(node.up.children[1].name)].endswith("_first"), \
                                "SA2 failure"
                            cur_range = int(node.up.children[1].name)

                        leaf_clade = SRangesClade(
                            frozenset({int(node.name)}),
                            ancestral_range=cur_range
                        )
                        clade_count_map[leaf_clade] += 1
                        observed_clades.add(leaf_clade)

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
        elif len(current_clade) == 0:
            # this is splitting off a range or SA
            continue

        # TODO new case that is not properly extracted:
        #  We have clade 84, 92 with anc 93 but there is no split recorded for it
        #  needs more fixing

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

            if current_split.descendant.clade and len(current_split.descendant.clade) == 1:
                leaf_observations = sum(
                    [clade_count_map[k] for k in leaf_clades if
                     k.clade == current_split.descendant.clade]
                )
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

    print("help me.")

    return None

# todo we will need to do a second post processing step over all ranges
#  that ended up being leaves without first and last, there add ranges only at the leaves
# todo maybe keep a dict of int --> is_range? or something like that?
