"""
Utilities for stratigraphic-range CCDs.

Provides data structures and algorithms for representing stratigraphic-range
clades, constructing stratigraphic-range CCDs from phylogenetic trees, and
CCD-MAP trees.
"""
import inspect
from collections import defaultdict
from dataclasses import dataclass

from brokilon.ccd.clades.sranges import SRangesClade
from brokilon.core import Tree

DEBUG_TREE = False


def get_callsite():
    """
    Get the source location of the caller's caller (for debugging not in use at the moment).

    Returns:
        A string containing the filename, line number, and function name.
    """
    frame = inspect.stack()[2]  # caller of caller
    return f"{frame.filename}:{frame.lineno} ({frame.function})"


@dataclass(frozen=True)
class AncestralSplit:
    """
    Represent an ancestral split in a stratigraphic-range CCD with orientation.

    Attributes:
        ancestor: The ancestral stratigraphic-range clade.
        descendant: The descendant stratigraphic-range clade.
    """
    ancestor: SRangesClade
    descendant: SRangesClade
    # source: str = field(default_factory=get_callsite)


def prelabel_tree(tree, taxon_map):
    """
    A function to prelabel a stratigraphic-range tree with differnt attributes.
    These attributes are node features called rangetype and sampledancestor to indicate whether the
    node is part of a range or displays a SA etc.

    :param tree: a stratigraphic-range tree
    :param taxon_map: the taxon map
    :return: opened_ranges, sampled_ancestors lists for both the ranges that were opened and the SAs
    """
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
            # non-regular leaf
            if current_taxon.endswith("_first"):
                other_leaves = [taxon_map[int(t.name)] for t in leaf.up if
                                t is not leaf]

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
                    leaf.up.add_feature("sampledancestor",
                                        taxon_map[int(leaf.name)][:-6])

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
    if DEBUG_TREE:
        print(
            f"Tree has {len(sampled_ancestors)} SAs and {len(opened_ranges)} ranges."
        )
    # returning set of all taxa that are a range and all that are SA
    return opened_ranges, sampled_ancestors


def get_sranges_map(trees, taxon_map):
    """Construct a stratigraphic-range CCD from phylogenetic trees.

    Processes the input trees to identify stratigraphic-range clades, ancestral
    splits, and sampled ancestors, and aggregates their counts to construct a
    stratigraphic-range CCD.

    :param trees: Collection of stratigraphic-range trees used to construct
        the CCD.
    :param taxon_map: Mapping between taxon labels and integer identifiers
        from the Nexus file.
    :return: A tuple containing the following four objects:

        - ``clade_count_map``: Dictionary mapping stratigraphic-range clades
          to their observed counts.
        - ``clade_split_count_map``: Nested dictionary mapping ancestral
          stratigraphic-range clades to their splits and observed counts.
        - ``sranges_set``: Set of stratigraphic ranges identified during
          processing.
        - ``sampled_ancestors``: Set of sampled ancestors identified in the
          trees.
    """
    clade_count_map = defaultdict(int)
    clade_split_count_map = defaultdict(lambda: defaultdict(int))

    # reverse_taxon_map = {value: key for key, value in taxon_map.items()}

    # for node in (node for t in trees for node in t.traverse("levelorder")):
    for t in trees:
        sranges_set, sampled_ancestors = prelabel_tree(t, taxon_map)

        for node in t.traverse("levelorder"):
            # We can ignore all leafs that are 0.0, i.e. internal leafs that encode a non leaf...
            if node.is_leaf():
                if node.dist == 0.0:
                    if getattr(node.up, "rangetype", None) == "leaf_range":
                        cur_range = None
                        if hasattr(node.up.up, "range"):
                            cur_range = f"{node.up.up.range}_first"
                        elif hasattr(node.up, "rangetype"):
                            if node.up.rangetype == "sampled_ancestor":
                                raise NotImplementedError("this might be necessary to do...")
                        current_leaf_range_start_clade = SRangesClade(
                            frozenset({f"{node.up.range}_first"}), cur_range)

                        clade_count_map[current_leaf_range_start_clade] += 1

                        # handle the empty clade as the descendant out of the range
                        # to indicate the start of the range
                        current_split = AncestralSplit(
                            descendant=SRangesClade(
                                frozenset({}),
                                f"{node.up.range}_first"
                            ),
                            ancestor=SRangesClade(
                                frozenset({f"{node.up.range}_first"}),
                                f"{node.up.range}_first"
                            )
                        )
                        clade_split_count_map[current_leaf_range_start_clade][
                            current_split] += 1
                    else:
                        continue
                if hasattr(node, "rangetype"):
                    if node.rangetype == "range_end":
                        if hasattr(node.up, "rangetype"):
                            if node.up.rangetype == "leaf_range":
                                # todo might have to add these clades too?
                                #  this would be to identify leaf ranges
                                continue
                        assert node.dist != 0.0, "If this happens we need to not allow it?!"
                        assert node.up.range == node.range, "The ranges don't match?"
                        assert (taxon_map[int(node.name)].replace("_last",
                                                                  "_first") ==
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
                    assert taxon_map[int(node.name)].endswith(
                        "_first"), "Failed leaf case..."
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
                        parent_clade_set = {
                            taxon_map[int(l.name)].replace("_last", "_first")
                            for l
                            in node}

                        child0_clade_set = {
                            taxon_map[int(l.name)].replace("_last", "_first")
                            for l
                            in node.children[0]}
                        child1_clade_set = {
                            taxon_map[int(l.name)].replace("_last", "_first")
                            for l
                            in node.children[1]}

                        assert parent_clade_set == child0_clade_set.union(
                            child1_clade_set)

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
                        clade_split_count_map[current_parent_clade][
                            current_split] += 1
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
                                        assert hasattr(node.up,
                                                       "sampledancestor"), \
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
                                {taxon_map[int(l.name)].replace("_last",
                                                                "_first") for l
                                 in node}
                            ),
                            current_range
                        )
                        clade_count_map[parent_clade] += 1
                    case "range_end":
                        # Keeping just the taxon name as range for now...
                        child0_clade_set = {
                            taxon_map[int(l.name)].replace("_last", "_first")
                            for l
                            in node.children[0]}
                        child1_clade_set = {
                            taxon_map[int(l.name)].replace("_last", "_first")
                            for l
                            in node.children[1]}
                        parent_clade_set = {
                            taxon_map[int(l.name)].replace("_last", "_first")
                            for l in node
                        }

                        assert parent_clade_set == child0_clade_set.union(
                            child1_clade_set), \
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
                            raise AssertionError(
                                f"Failure: range end should have "
                                f"single taxon leaf set which is the range...")
                        parent_range = None
                        parent_clade_should_be_counted = True
                        if node.up:
                            if hasattr(node.up, "rangetype"):
                                match node.up.rangetype:
                                    # todo this case with the the case elif below? redundancy?
                                    case "range_end":
                                        parent_range = f"{node.up.range}_first"
                                    case "sampled_ancestor":
                                        assert hasattr(node.up,
                                                       "sampledancestor"), \
                                            "This needs fixing"
                                        parent_range = f"{node.up.sampledancestor}_first"
                                    case "range_start":
                                        # We don't need to count this parent
                                        # it was counted already at the start of this range
                                        # and nothing happend inside the range
                                        parent_clade_should_be_counted = False
                                        assert node.up.range == node.range, \
                                            "something wrong with node.up.up case"
                                        if hasattr(node.up.up, "range"):
                                            parent_range = f"{node.up.up.range}_first"
                                        elif hasattr(node.up.up,
                                                     "sampledancestor"):
                                            parent_range = f"{node.up.up.sampledancestor}_first"
                                    case _:
                                        pass
                            elif hasattr(node.up, "range"):
                                parent_range = f"{node.up.range}_first"

                        parent_clade = SRangesClade(
                            frozenset(parent_clade_set),
                            parent_range
                        )
                        # here we are closing a range so we need to count the split
                        # if the parent is also the start of the range we don't count it
                        # otherwise the parent is not also the start of the range, i.e.
                        # something else happens inside the range (range-birth) and the parent clade
                        # needs to be counted
                        if parent_clade_should_be_counted:
                            clade_count_map[parent_clade] += 1

                        clade_split_count_map[parent_clade][current_split] += 1
                    case "leaf_range":
                        # todo this seems to be counting leaf ranges twice with the new version
                        #  hence it will likely be removed...

                        # We are looking at a node that is parent of a range that is a leaf
                        # We don't need to add a split here?
                        # We only count the clade, which is the taxon and itself as a range?

                        # # todo fix this to something less messy...
                        # current_leaf = [l.name for l in node if l.dist == 0.0][
                        #     0]
                        # range_taxon = reverse_taxon_map[f"{node.range}_first"]
                        #
                        # # todo what is this assert doing?
                        # assert range_taxon == int(current_leaf), \
                        #     (f"Failure: {range_taxon} != {current_leaf}")
                        #
                        # # Leaf ranges can not have themselves as a range,
                        # # but there can be other ranges as ancestors above
                        # # Leaf ranges have no ancestral range of themselves because
                        # # otherwise this information is assumed to be present above the leaf
                        # # and therefore breaks the logic
                        # # One way to avoid this and fully encode leaf ranges is
                        # # to have _first and _last taxon as part of the clades...
                        # current_range = None
                        # if hasattr(node.up, "range"):
                        #     current_range = f"{node.up.range}_first"
                        # elif hasattr(node.up, "sampledancestor"):
                        #     current_range = f"{node.up.sampledancestor}_first"
                        # clade_count_map[
                        #     SRangesClade(
                        #         frozenset({taxon_map[int(current_leaf)]}),
                        #         current_range
                        #     )
                        # ] += 1
                        continue
                    case _:
                        raise ValueError("Unsupported value...")

            else:
                parent_clade_set = {
                    taxon_map[int(l.name)].replace("_last", "_first") for l in
                    node
                }

                child0_clade_set = {
                    taxon_map[int(l.name)].replace("_last", "_first") for l in
                    node.children[0]
                }
                child1_clade_set = {
                    taxon_map[int(l.name)].replace("_last", "_first") for l in
                    node.children[1]
                }

                assert child1_clade_set.union(
                    child0_clade_set) == parent_clade_set, \
                    "Something is wrong with the clades..."
                cur_range = None
                parent_range = None
                range_birth = False
                if hasattr(node, "range"):
                    range_birth = True
                    cur_range = f"{node.range}_first"
                if node.up:
                    if hasattr(node.up, "rangetype"):
                        match node.up.rangetype:
                            case "range_end":
                                parent_range = f"{node.up.range}_first"
                            case "sampled_ancestor":
                                assert hasattr(node.up,
                                               "sampledancestor"), "This needs fixing"
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

                parent_clade = SRangesClade(frozenset(parent_clade_set),
                                            parent_range)

                if node.children[0].orientation == "ancestor":
                    current_split = AncestralSplit(
                        ancestor=SRangesClade(frozenset(child0_clade_set),
                                              cur_range),
                        descendant=SRangesClade(frozenset(child1_clade_set),
                                                cur_range),
                    )
                else:
                    current_split = AncestralSplit(
                        ancestor=SRangesClade(frozenset(child1_clade_set),
                                              cur_range),
                        descendant=SRangesClade(frozenset(child0_clade_set),
                                                cur_range),
                    )
                if range_birth and getattr(node.up, "rangetype", "nothing") == "range_start":
                    # this parent clade was already added with the start range and does not
                    # need to be counted again (double count issue)
                    if DEBUG_TREE:
                        print(f"[DEBUG]: First range birth in range, no clade count added..")
                else:
                    clade_count_map[parent_clade] += 1
                clade_split_count_map[parent_clade][current_split] += 1

    return (
        dict(clade_count_map),
        {clade: dict(splits) for clade, splits in clade_split_count_map.items()},
        sranges_set, sampled_ancestors
    )


def get_sranges_map_tree(
        clade_count_map,
        clade_split_count_map,
        sranges_set, sampled_ancestors
):
    """Construct the stratigraphic-range CCD-MAP tree.

    Resolves clades in order of increasing size using their observed counts
    and ancestral split counts. The resolved clades are then used in a
    dynamic programming algorithm to construct the MAP tree, incorporating
    stratigraphic ranges and sampled ancestors.

    Args:
        clade_count_map: Dictionary mapping stratigraphic-range clades to
            their observed counts.
        clade_split_count_map: Nested dictionary mapping ancestral clades to
            their ancestral splits and corresponding counts.
        sranges_set: Set of stratigraphic ranges identified during CCD
            construction.
        sampled_ancestors: Set of sampled ancestors identified during CCD
            construction.

    Returns:
        A tuple containing:

        - `seen_resolved_clades`: Dictionary of clades resolved during MAP
          tree construction, primarily for debugging.
        - `map_tree`: The resulting stratigraphic-range CCD-MAP tree.
    """
    seen_resolved_clades = {}
    for current_clade in sorted(clade_count_map.keys(), key=len):
        if len(current_clade) < 2:
            # leaf clades or SA/range specific clades to skip
            continue

        for current_split in clade_split_count_map[current_clade]:

            # anc_prob, desc_prob = 0, 0

            if len(current_split.ancestor.clade) == 1:
                anc_prob = 1
            elif len(current_split.ancestor.clade) == 0:
                anc_prob = 1
            else:
                # [c for c in seen_resolved_clades if c.clade == current_split.ancestor.clade],
                # current_split.ancestor, seen_resolved_clades[current_split.ancestor]
                anc_prob = seen_resolved_clades[current_split.ancestor][0]

            if len(current_split.descendant.clade) == 1:
                desc_prob = 1
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
                    if DEBUG_TREE:
                        print(
                            f"DEBUG[SRB]: {current_clade} is being set to {current_split} "
                            f"with probability {split_probability}"
                        )
                    seen_resolved_clades[current_clade] = (split_probability, current_split, False)
                elif seen_resolved_clades[current_clade][0] == split_probability:
                    # Tie breaking randomly, either keep the old or pick the new
                    import random
                    if random.random() < 0.5:
                        chosen_prob, chosen_split = split_probability, current_split
                    else:
                        chosen_prob, chosen_split = seen_resolved_clades[
                            current_clade][:2]

                    if DEBUG_TREE:
                        print(
                            f"DEBUG[SRB]: {current_clade} is being set to {current_split} "
                            f"with probability {split_probability}"
                        )
                    # resolving tiebreak in the seen_resolved_clades, True indicates the tiebreak
                    seen_resolved_clades[current_clade] = (chosen_prob, chosen_split, True)
            else:
                if DEBUG_TREE:
                    print(
                        f"DEBUG[SRB]: {current_clade} is being set to {current_split} "
                        f"with probability {split_probability}"
                    )
                seen_resolved_clades[current_clade] = (split_probability, current_split, False)
    # End of seen_resolved_clades construction

    # output = {}
    max_key_value = len(max(seen_resolved_clades.keys()))
    all_root_clades = [k for k in seen_resolved_clades if
                       len(k) == max_key_value]
    max_root_prob = max(
        seen_resolved_clades[root][0] for root in all_root_clades)
    best_roots = [root for root in all_root_clades
                  if seen_resolved_clades[root][0] == max_root_prob]

    if len(best_roots) > 1:
        raise NotImplementedError(
            "Currently no tie breaking support for multiple roots")
    else:
        if DEBUG_TREE:
            print(f"DEBUG[root]: {seen_resolved_clades[best_roots[0]]}")
        root = best_roots[0]

    map_tree = get_sranges_tree_from_seen_resolved_clades(
        seen_resolved_clades,
        root,
        sranges_set, sampled_ancestors
    )
    #
    return seen_resolved_clades, map_tree


def get_sranges_tree_from_seen_resolved_clades(
        seen_resolved_clades,
        root_clade,
        sranges_set, sampled_ancestors
):
    """Construct a stratigraphic-range tree from resolved clades.

    Recursively builds a tree from the clades and ancestral splits resolved
    during MAP tree construction. Adds internal nodes, stratigraphic-range
    nodes, and sampled ancestors where appropriate, preserving their
    ancestral-descendant relationships and orientation annotations.

    :param seen_resolved_clades: Mapping of clades and their resolved ancestral
        splits used to construct the tree.
    :param root_clade: Root clade of the stratigraphic-range CCD-MAP tree.
    :param sranges_set: Set of stratigraphic ranges identified during CCD
        construction.
    :param sampled_ancestors: Set of sampled ancestors identified during CCD
        construction.
    :return: The constructed tree, including its stratigraphic-range and
        orientation annotations.
    """
    # todo ignoring tie breaking for now...

    out_tree = Tree(
        dist=1.0,
        name="root"
    )
    if root_clade.ancestral_range in sranges_set:
        out_tree.add_feature("ancestral_range", root_clade.ancestral_range)
    out_tree.add_feature("orientation", "ancestor")

    icount = 1

    def recursive_sranges_children(node, split):
        """
        Recursively add the descendant clades of an ancestral split.

        Processes the ancestral and descendant clades of a split and adds their
        corresponding nodes and relationships to the output tree.
        """
        nonlocal icount, seen_resolved_clades
        nonlocal sranges_set

        if DEBUG_TREE:
            print(f"DEBUG: {split} is being added")

        cur_anc = split.ancestor
        cur_descendant = split.descendant
        if split.descendant.clade == frozenset():
            # todo this should not be fixed here, but might help me for now
            cur_anc = split.descendant
            cur_descendant = split.ancestor

        def add_clade(parent_node, clade, orientation, second_part_of_split):
            """
            Add a clade and its associated nodes to the output tree.

            Adds the clade beneath the specified parent node, handling its orientation,
            position within the split, and associated stratigraphic-range nodes.
            """
            nonlocal icount, seen_resolved_clades
            nonlocal sranges_set

            def find_range_children(parent):
                """
                Find the previously constructed internal/range-node pair
                below one of the parent's children.
                """
                for child in parent.children:
                    try:
                        adjusted_parent = next(
                            c for c in child.children
                            if 'internal_' in c.name
                        )
                        start_range = next(
                            c for c in child.children
                            if 'internal_' not in c.name
                        )
                        return adjusted_parent, start_range

                    except StopIteration:
                        continue

                raise NotImplementedError(
                    "Could not find expected internal/range children."
                )

            # ------------------------------------------------------------------
            # Empty/special clade
            # ------------------------------------------------------------------
            if not clade.clade:
                # Special clade to encode ranges and SAs
                parent_range = getattr(parent_node, 'ancestral_range', None)

                assert clade.ancestral_range not in second_part_of_split.clade, (
                    "If not true, we are not closing a range."
                )

                match parent_range:
                    case clade.ancestral_range:
                        # Parent range is same as current node:
                        # we are adding a range end node.
                        range_end_node = parent_node.add_child(
                            name=clade.ancestral_range.replace("_first", "_last"),
                            dist=1.662 if DEBUG_TREE else 1
                        )
                        range_end_node.add_feature(
                            'orientation',
                            orientation
                        )
                        range_end_node.add_feature(
                            'ancestral_range',
                            clade.ancestral_range
                        )

            # ------------------------------------------------------------------
            # Leaf
            # ------------------------------------------------------------------
            elif len(clade.clade) == 1:
                label = next(iter(clade.clade))

                if (
                        clade.ancestral_range not in clade.clade
                        and second_part_of_split.clade != frozenset()
                ):
                    # Adding a leaf after the end of a range,
                    # no need to add internal special nodes.
                    if (
                            clade.ancestral_range
                            and len(parent_node.children) == 2
                    ):
                        # Need to adjust the leaf parent node.
                        adjusted_parent_node, start_range_node = (
                            find_range_children(parent_node)
                        )
                        if label in sranges_set:
                            # have to add a leaf range

                            range_start_parent = adjusted_parent_node.add_child(
                                name=f"internal_leafrange",
                                dist=7.501 if DEBUG_TREE else 1,
                            )
                            range_start_parent.add_feature("orientation", orientation)

                            range_start_leaf = range_start_parent.add_child(
                                name=label,
                                dist=0.157 if DEBUG_TREE else 0,
                            )
                            range_start_leaf.add_feature("orientation", "descendant")

                            leaf = range_start_parent.add_child(
                                name=label.replace("_first", "_last"),
                                dist=1.145 if DEBUG_TREE else 1,
                            )
                            leaf.add_feature("orientation", "ancestor")
                            leaf.add_feature("ancestral_range", label)

                        else:
                            leaf = adjusted_parent_node.add_child(
                                name=label,
                                dist=1.718 if DEBUG_TREE else 1,
                            )
                            leaf.add_feature("orientation", orientation)
                    elif (
                            clade.ancestral_range in second_part_of_split.clade
                            and getattr(parent_node, "ancestral_range",
                                        None) != clade.ancestral_range
                    ):
                        adjusted_parent_node, start_range_node = find_range_children(parent_node)

                        if label in sranges_set:
                            # have to add a leaf range

                            range_start_parent = adjusted_parent_node.add_child(
                                name=f"internal_leafrange",
                                dist=1.733 if DEBUG_TREE else 1,
                            )
                            range_start_parent.add_feature("orientation", orientation)

                            range_start_leaf = range_start_parent.add_child(
                                name=label,
                                dist=0.145 if DEBUG_TREE else 0,
                            )
                            range_start_leaf.add_feature("orientation", "descendant")

                            leaf = range_start_parent.add_child(
                                name=label.replace("_first", "_last"),
                                dist=1.34 if DEBUG_TREE else 1,
                            )
                            leaf.add_feature("orientation", "ancestor")
                            leaf.add_feature("ancestral_range", label)

                        else:

                            leaf = adjusted_parent_node.add_child(
                                name=label,
                                dist=1.754 if DEBUG_TREE else 1,
                            )
                            leaf.add_feature("orientation", orientation)
                    else:

                        if label in sranges_set:
                            # have to add a leaf range

                            range_start_parent = parent_node.add_child(
                                name=f"internal_leafrange",
                                dist=1.764 if DEBUG_TREE else 1,
                            )
                            range_start_parent.add_feature("orientation", orientation)

                            range_start_leaf = range_start_parent.add_child(
                                name=label,
                                dist=0.124 if DEBUG_TREE else 0,
                            )
                            range_start_leaf.add_feature("orientation", "descendant")

                            leaf = range_start_parent.add_child(
                                name=label.replace("_first", "_last"),
                                dist=1.12 if DEBUG_TREE else 1,
                            )
                            leaf.add_feature("orientation", "ancestor")
                            leaf.add_feature("ancestral_range", label)

                        else:
                            leaf = parent_node.add_child(
                                name=label,
                                dist=1.784 if DEBUG_TREE else 1,
                            )
                            leaf.add_feature("orientation", "descendant")

                elif (
                        clade.ancestral_range
                        and getattr(second_part_of_split, 'ancestral_range', None)
                        == getattr(clade, 'ancestral_range', None)
                        and getattr(parent_node, 'ancestral_range', None)
                        != getattr(clade, 'ancestral_range', None)
                ):
                    if clade.ancestral_range not in clade.clade:
                        if clade.ancestral_range in sampled_ancestors:
                            # this case we are adding a leaf that has a range on top
                            sa_parent = parent_node.add_child(
                                name=f"internal_{icount}",
                                dist=1.8 if DEBUG_TREE else 1,
                            )
                            sa_parent.add_feature(
                                "orientation",
                                orientation
                            )
                            sa_leaf = sa_parent.add_child(
                                name=clade.ancestral_range,
                                dist=1.808 if DEBUG_TREE else 0,
                            )
                            sa_leaf.add_feature('orientation', "ancestor")

                            if label in sranges_set:
                                # have to add a leaf range

                                range_start_parent = sa_parent.add_child(
                                    name=f"internal_leafrange",
                                    dist=1.90 if DEBUG_TREE else 1,
                                )
                                range_start_parent.add_feature("orientation", orientation)

                                range_start_leaf = range_start_parent.add_child(
                                    name=label,
                                    dist=0.1354 if DEBUG_TREE else 0,
                                )
                                range_start_leaf.add_feature("orientation", "descendant")

                                leaf = range_start_parent.add_child(
                                    name=label.replace("_first", "_last"),
                                    dist=1.986 if DEBUG_TREE else 1,
                                )
                                leaf.add_feature("orientation", "ancestor")
                                leaf.add_feature("ancestral_range", label)

                            else:
                                leaf = sa_parent.add_child(
                                    name=label,
                                    dist=1.987 if DEBUG_TREE else 1,
                                )
                                leaf.add_feature("orientation", "descendant")

                        else:
                            # todo fix this, why is the assert now wrong, is it because the bugfix?
                            assert clade.ancestral_range in sranges_set, "should be the case"
                            # this case we are adding a leaf that has a range on top
                            start_range_parent = parent_node.add_child(
                                name=f"internal_{icount}",
                                dist=10 if DEBUG_TREE else 1,
                            )
                            start_range_parent.add_feature("orientation", orientation)

                            range_start_leaf = start_range_parent.add_child(
                                name=clade.ancestral_range,
                                dist=2.4 if DEBUG_TREE else 0,
                            )
                            range_start_leaf.add_feature('orientation', "descendant")

                            extra_internal_node = start_range_parent.add_child(
                                name=f"internal_extra_{icount}",
                                dist=3 if DEBUG_TREE else 1,
                            )
                            # if getattr(clade, 'ancestral_range', None) in clade.clade:
                            extra_internal_node.add_feature('ancestral_range',
                                                            clade.ancestral_range)
                            extra_internal_node.add_feature('orientation', "ancestor")

                            # assuming the range here is always ontop of the leaf,
                            # never actually a leaf itself, therefore dist 0.0
                            range_end = extra_internal_node.add_child(
                                name=clade.ancestral_range.replace('_first', '_last'),
                                dist=0.868 if DEBUG_TREE else 0,
                            )
                            range_end.add_feature("orientation", "ancestor")
                            range_end.add_feature("ancestral_range", clade.ancestral_range)

                            if label in sranges_set:
                                # have to add a leaf range

                                range_start_parent = extra_internal_node.add_child(
                                    name=f"internal_leafrange",
                                    dist=1.90 if DEBUG_TREE else 1,
                                )
                                range_start_parent.add_feature("orientation", orientation)

                                range_start_leaf = range_start_parent.add_child(
                                    name=label,
                                    dist=0.1354 if DEBUG_TREE else 0,
                                )
                                range_start_leaf.add_feature("orientation", "descendant")

                                leaf = range_start_parent.add_child(
                                    name=label.replace("_first", "_last"),
                                    dist=1.986 if DEBUG_TREE else 1,
                                )
                                leaf.add_feature("orientation", "ancestor")
                                leaf.add_feature("ancestral_range", label)

                            else:
                                leaf = extra_internal_node.add_child(
                                    name=label,
                                    dist=1.9823 if DEBUG_TREE else 1,
                                )
                                leaf.add_feature("orientation", "descendant")

                    elif label in sranges_set:
                        assert clade.ancestral_range in clade.clade, "would hope so..."
                        start_range_parent = parent_node.add_child(
                            name=f"internal_{icount}",
                            dist=10 if DEBUG_TREE else 1,
                        )
                        start_range_parent.add_feature("orientation", orientation)

                        range_start_leaf = start_range_parent.add_child(
                            name=label,
                            dist=2.4 if DEBUG_TREE else 0,
                        )
                        range_start_leaf.add_feature(
                            'orientation',
                            "descendant"
                        )

                        extra_internal_node = start_range_parent.add_child(
                            name=f"internal_extra_{icount}",
                            dist=3 if DEBUG_TREE else 1,
                        )
                        if getattr(clade, 'ancestral_range', None) in clade.clade:
                            extra_internal_node.add_feature(
                                'ancestral_range',
                                clade.ancestral_range
                            )
                        extra_internal_node.add_feature(
                            'orientation',
                            "ancestor"
                        )

                        leaf = extra_internal_node.add_child(
                            name=label.replace('_first', '_last'),
                            dist=1.997 if DEBUG_TREE else 1,
                        )
                        leaf.add_feature("orientation", "ancestor")
                        leaf.add_feature("ancestral_range", clade.ancestral_range)
                    else:
                        raise NotImplementedError(f"sa case with {label}")
                else:
                    if label in sranges_set and clade.ancestral_range != label:

                        start_range_parent = parent_node.add_child(
                            name=f"internal_{icount}",
                            dist=1.023 if DEBUG_TREE else 1,
                        )
                        start_range_parent.add_feature("orientation", "descendant")

                        range_start_leaf = start_range_parent.add_child(
                            name=label,
                            dist=0.45 if DEBUG_TREE else 0,
                        )
                        range_start_leaf.add_feature('orientation', "descendant")

                        leaf_range_end = start_range_parent.add_child(
                            name=label.replace('_first', '_last'),
                            dist=4.24 if DEBUG_TREE else 0,
                        )
                        leaf_range_end.add_feature("orientation", "ancestor")
                        leaf_range_end.add_feature("ancestral_range", clade.ancestral_range)

                    else:
                        label = label.replace('_first', '_last')

                        leaf = parent_node.add_child(
                            name=label,
                            dist=10 if DEBUG_TREE else 1,
                        )

                        leaf.add_feature("orientation", orientation)
                        if next(iter(clade.clade)) == clade.ancestral_range:
                            leaf.add_feature('ancestral_range', clade.ancestral_range)

            # ------------------------------------------------------------------
            # Internal node
            # ------------------------------------------------------------------
            else:
                if (
                        clade.ancestral_range and
                        clade.ancestral_range != getattr(parent_node, 'ancestral_range', None)
                ):
                    if clade.ancestral_range in clade.clade:
                        start_range_parent = parent_node.add_child(
                            name=f"internal_{icount}",
                            dist=6 if DEBUG_TREE else 1,
                        )
                        start_range_parent.add_feature(
                            'orientation',
                            orientation
                        )

                        start_range_leaf = start_range_parent.add_child(
                            name=clade.ancestral_range,
                            dist=0.123 if DEBUG_TREE else 0,
                        )
                        start_range_leaf.add_feature('orientation', "ancestor")

                        extra_internal_node = start_range_parent.add_child(
                            name=f"internal_extra_{icount}",
                            dist=3 if DEBUG_TREE else 1,
                        )
                        if hasattr(clade, 'ancestral_range'):
                            extra_internal_node.add_feature(
                                'ancestral_range',
                                clade.ancestral_range
                            )
                        extra_internal_node.add_feature(
                            'orientation',
                            "ancestor"
                        )

                        internal_node = extra_internal_node.add_child(
                            name=f"internal_{icount}",
                            dist=23 if DEBUG_TREE else 1,
                        )
                        if clade.ancestral_range in clade.clade:
                            internal_node.add_feature(
                                'ancestral_range',
                                clade.ancestral_range
                            )
                        internal_node.add_feature(
                            "orientation",
                            "descendant"
                        )

                    else:
                        if (
                                clade.ancestral_range in sranges_set
                                and clade.ancestral_range
                                not in second_part_of_split.clade
                        ):
                            range_start_leaf = parent_node.add_child(
                                name=clade.ancestral_range,
                                dist=2.443 if DEBUG_TREE else 0,
                            )
                            range_start_leaf.add_feature(
                                'orientation',
                                "descendant"
                            )

                            internal_extra_node = parent_node.add_child(
                                name=f"internal_extra_{icount}",
                                dist=2.443 if DEBUG_TREE else 1
                            )
                            if hasattr(clade, 'ancestral_range'):
                                internal_extra_node.add_feature(
                                    'ancestral_range',
                                    clade.ancestral_range
                                )
                            internal_extra_node.add_feature(
                                'orientation',
                                "ancestor"
                            )

                            range_end_node = internal_extra_node.add_child(
                                name=clade.ancestral_range.replace(
                                    "_first",
                                    "_last"
                                ),
                                dist=2.443 if DEBUG_TREE else 0,
                            )
                            range_end_node.add_feature(
                                'orientation',
                                "ancestor"
                            )
                            if hasattr(clade, 'ancestral_range'):
                                range_end_node.ancestral_range = (
                                    clade.ancestral_range
                                )

                            internal_node = internal_extra_node.add_child(
                                name=f"internal_{icount}",
                                dist=1.4675 if DEBUG_TREE else 1
                            )
                            internal_node.add_feature("orientation", "descendant")
                        else:
                            if clade.ancestral_range in sranges_set:

                                adjusted_parent_node, start_range_node = (
                                    find_range_children(parent_node)
                                )

                                internal_node = (
                                    adjusted_parent_node.add_child(
                                        name=f"internal_{icount}",
                                        dist=4 if DEBUG_TREE else 1,
                                    )
                                )

                                internal_node.add_feature(
                                    'orientation',
                                    "descendant"
                                )

                            else:
                                # SA case is being added.
                                sampled_ancestor_node = parent_node.add_child(
                                    name=clade.ancestral_range,
                                    dist=0.144 if DEBUG_TREE else 0,
                                )
                                sampled_ancestor_node.add_feature(
                                    'orientation',
                                    orientation
                                )

                                internal_node = parent_node.add_child(
                                    name=f"internal_{icount}",
                                    dist=4.5 if DEBUG_TREE else 1,
                                )

                                internal_node.add_feature('orientation', "descendant")

                else:
                    internal_node = parent_node.add_child(
                        name=f"internal_{icount}",
                        dist=7 if DEBUG_TREE else 1,
                    )

                    internal_node.add_feature("orientation", orientation)
                    if getattr(clade, "ancestral_range", None) in clade.clade:
                        internal_node.add_feature(
                            "ancestral_range",
                            clade.ancestral_range
                        )

                icount += 1

                # Recurse if this clade was resolved.
                if clade in seen_resolved_clades:

                    if DEBUG_TREE:
                        for s in seen_resolved_clades[clade]:
                            print(f"DEBUG[all splits]: {s}")

                    child_split = seen_resolved_clades[clade][1]
                    recursive_sranges_children(
                        internal_node,
                        child_split
                    )

        add_clade(node, cur_anc, "ancestor", cur_descendant)
        add_clade(node, cur_descendant, "descendant", cur_anc)

    recursive_sranges_children(
        out_tree,
        seen_resolved_clades[root_clade][1]
    )

    return out_tree
