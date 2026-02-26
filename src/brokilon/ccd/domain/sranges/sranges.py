from collections import defaultdict

import logging
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
            ranges_clade.add(int(leaf.name))
    return ranges_clade


def get_sranges_map(trees, taxon_map, ccd_type=1):

    clade_count_map = defaultdict(int)
    clade_split_count_map = defaultdict(lambda: defaultdict(int))

    reverse_taxon_map = {value: key for key, value in taxon_map.items()}

    for node in (node for t in trees for node in t.traverse("levelorder")):
        if len(node.children) == 1:
            logging.debug(f"Found a node with only 1 child")
            continue

        if len(node) > 1:

            # Encoding ranged, int if starting a range with child 1 or 0
            # It will be "End" if we are ending a range with this node
            # None for all the inbetween, i.e. inside and outside ranges

            current_range = None

            # todo with the new ranges_clade function i think I am adding ranges twice in some cases.

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
                        # {int(leaf.name) for leaf in
                        #           if
                        #          leaf.dist != 0.0}
                    # c1_leaves = {int(leaf.name) for leaf in node.children[1] if
                    #              leaf.dist != 0.0}

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
                        ancestor= anc_clade,
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

                # todo this is wrong, need to have ranges in there too, only once though..
                leaves = get_ranges_clade(node.children[current_range], taxon_map)
                    # {int(leaf.name) for leaf in node.children[current_range] if
                    #          leaf.dist != 0.0}

                # adding the range using the first occurence of the fossil
                leaves.add(reverse_taxon_map[f"{node.children[current_range].range}_first"])
                current_range_start_clade = SRangesClade(
                    frozenset(leaves),
                    ancestral_range=None
                )
                clade_count_map[current_range_start_clade] += 1
            else:
                raise ValueError(f"Unexpected value for current range: {current_range}!")

        elif node.is_leaf():
            if node.dist == 0.0:
                # todo this might cause problems if we ignore leaves that are ending a range?...
                # todo could check with node.range_last == node.name that would be such a case?
                # we can ignore these because they are handled above
                continue
            leaf_ancestor_range = None
            if hasattr(node.up, "range"):
                leaf_ancestor_range = reverse_taxon_map[f"{node.up.range}_first"]
            clade_count_map[SRangesClade(
                frozenset({int(node.name)}), ancestral_range=leaf_ancestor_range
            )] += 1

    # todo do not return default dict, otherwise its easy to add things etc...
    return clade_count_map, clade_split_count_map