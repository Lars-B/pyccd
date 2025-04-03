import os.path
import re
import warnings
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
from docutils.io import InputError

from pyccd.read_breath_nexus import read_transmission_nexus_history


@dataclass(frozen=True)
class TransmissionClade:
    """Class to keep clades with additional transmission block information"""
    clade: frozenset
    has_block: bool
    # blockcount: int  # not sure if I want this in here or not...

    def __len__(self):
        return len(self.clade)

    def __lt__(self, other):
        if isinstance(other, TransmissionClade):
            return len(self.clade) < len(other.clade)
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, TransmissionClade):
            return len(self.clade) > len(other.clade)
        return NotImplemented


def get_transmission_clades(tree, blockcountmap):
    treestr = tree.write(format=8)
    clades = set()
    if not (treestr[0] == '(' and treestr[-2] == ')' and treestr[-1] == ';'):
        # todo the above second statement is not true: there might be root info that we can ignore
        raise InputError("Invalid tree string given! (no ';' at the end)")

    # Adding leafs separately
    leafnamepattern = r"%(-?\d+/-?\d+)%"
    # the above regex assumes that the taxa are integers and internal nodes do not have an
    # interger in the first part of the block.
    matches = re.findall(leafnamepattern, treestr)
    for m in matches:
        taxa, blockcount = map(int, m.split("/"))
        if blockcount == -1:
            leafclade = TransmissionClade(frozenset({taxa}), False)
        else:
            leafclade = TransmissionClade(frozenset({taxa}), True)
            if leafclade in blockcountmap:
                blockcountmap[leafclade] += blockcount
            else:
                blockcountmap[leafclade] = blockcount
        clades.add(leafclade)

    leafreplacementpattern = r"%(-?\d+)/-?\d+%"

    def replace_match(match):
        return match.group(1)  # Extract leafname from above pattern (the part before the slash)

    # Replace all matches with integer1
    treestr = re.sub(leafreplacementpattern, replace_match, treestr)

    # Adding all non leaf clades (except the root)
    opend = []
    re_brackets_internals = re.compile(r"\(|\)|%/-?\d+%")
    internalnodenamepattern = re.compile(r"%/-?\d+%")
    for i in range(1, len(treestr) - 2):
        if treestr[i] == '(':
            opend.append(i)
        elif treestr[i] == ')':
            if not opend:
                raise InputError("Invalid tree string given! (to many ')')")
            # current_string = treestr[opend[-1]:i]

            match_interenalnodename = re.search(internalnodenamepattern, treestr[i+1:])
            if not match_interenalnodename:
                raise ValueError("Something went wrong looking for the next internal nodes name")
            blockcount = int(
                match_interenalnodename.group(0).replace("%", "").split("/")[1]
            )
            current_clade_set = frozenset(
                re.sub(re_brackets_internals, "", treestr[opend[-1]:i]).split(",")
            )
            if blockcount == -1:
                curclade = TransmissionClade(current_clade_set, False)
            else:
                curclade = TransmissionClade(current_clade_set, True)
                if curclade in blockcountmap:
                    blockcountmap[curclade].append(blockcount)
                else:
                    blockcountmap[curclade] = [blockcount]
            clades.add(curclade)
            del opend[-1]
    if opend:
        raise InputError("Invalid tree string given! (to many '(')")
    return clades


def get_transmission_maps(trees):

    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled

    blockcount_map = {}
    branch_lengths_map = {}

    for _, t in enumerate(trees):
        for node in t.traverse("levelorder"):
            if len(node) > 1:
                c = node.children
                c0_leafs = set()
                for leaf in c[0]:
                    # c0_leafs.add(int(leaf.name.replace("%", "").split("/")[0]))
                    c0_leafs.add(int(leaf.name))
                c1_leafs = set()
                for leaf in c[1]:
                    # c1_leafs.add(int(leaf.name.replace("%", "").split("/")[0]))
                    c1_leafs.add(int(leaf.name))

                parent_clade_set = frozenset(sorted(c0_leafs.union(c1_leafs)))

                if node.name:
                    # blockcount = int(node.name.replace("%", "").split("/")[1])
                    assert hasattr(node, "blockcount"), ("The nodes should have the blockcount "
                                                         "attribute!")
                    blockcount = node.blockcount
                    if blockcount == -1:
                        parent_clade = TransmissionClade(parent_clade_set, False)
                    else:
                        parent_clade = TransmissionClade(parent_clade_set, True)
                        if parent_clade in blockcount_map:
                            blockcount_map[parent_clade].append(blockcount)
                        else:
                            blockcount_map[parent_clade] = [blockcount]
                else:
                    # this is the root, apparently we don't have the name atm
                    parent_clade = TransmissionClade(parent_clade_set, False)

                # adding distance of parent clade to map of branch lengths
                if parent_clade in branch_lengths_map:
                    branch_lengths_map[parent_clade].append(node.dist)
                else:
                    branch_lengths_map[parent_clade] = [node.dist]

                m1[parent_clade] += 1

                blockcount_c0 = c[0].blockcount
                blockcount_c1 = c[1].blockcount
                child0_clade = TransmissionClade(frozenset(c0_leafs), blockcount_c0 != -1)
                child1_clade = TransmissionClade(frozenset(c1_leafs), blockcount_c1 != -1)

                # in m2 the split clade with the lower int for taxa is entered first
                if min(c0_leafs) < min(c1_leafs):
                    m2[(parent_clade, child0_clade, child1_clade)] += 1
                else:
                    m2[(parent_clade, child1_clade, child0_clade)] += 1
            elif len(node) == 1:
                # leaf node for which we need to add the blockcount to the blockcount_map
                leaf_label = int(node.name)
                blockcount = node.blockcount

                leaf_clade = TransmissionClade(frozenset({leaf_label}), (blockcount != -1))
                # adding branch length of the leaf node to the dict
                if leaf_clade in branch_lengths_map:
                    branch_lengths_map[leaf_clade].append(node.dist)
                else:
                    branch_lengths_map[leaf_clade] = [node.dist]

                # if the clade has a block we need to keep track of the blockcount for summaries
                if blockcount != -1:
                    if leaf_clade in blockcount_map:
                        blockcount_map[leaf_clade].append(blockcount)
                    else:
                        blockcount_map[leaf_clade] = [blockcount]
    return m1, m2, blockcount_map, branch_lengths_map


def get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map):
    seen_resolved_clades = {}
    # sorted list of all clades, small (cherries) to big
    all_clades = sorted(list(m1.keys()), key=len)

    for current_clade in all_clades:
        # the following are all triplets that represent how the current clade splits
        all_splits_current = [i for i in m2 if i[0] == current_clade]

        for current_split in all_splits_current:
            child1, child2 = current_split[1], current_split[2]
            # c1_prob = 0
            # c2_prob = 0

            if len(child1) == 1:
                # it is a leaf, hence probability is 1
                c1_prob = 1
            else:
                if child1 in seen_resolved_clades:
                    c1_prob = seen_resolved_clades[child1][0]
                else:
                    raise ValueError("Should never get here?!")

            if len(child2) == 1:
                # it is a leaf, hence probability is 1
                c2_prob = 1
            else:
                if child2 in seen_resolved_clades:
                    c2_prob = seen_resolved_clades[child2][0]
                else:
                    raise ValueError("Should never get here?!")

            cur_prob = m2[current_split] / m1[current_split[0]]
            split_prob = c1_prob * c2_prob * cur_prob

            if current_split[0] in seen_resolved_clades:
                if seen_resolved_clades[current_split[0]][0] < split_prob:
                    seen_resolved_clades[current_split[0]] = (split_prob, current_split)
                elif seen_resolved_clades[current_split[0]][0] == split_prob:
                    # todo add a way to make the tie breaking deterministic or seedable
                    warnings.warn("Tie breaking in effect.")
                    seen_resolved_clades[current_split[0]] = (split_prob, current_split)
            else:
                seen_resolved_clades[current_split[0]] = (split_prob, current_split)

    root_clade = max(seen_resolved_clades.keys())
    working_list = [root_clade]
    output = {}

    while working_list:
        cur_parent = working_list.pop()
        # output.append(seen_resolved_clades[cur_parent][1])
        cur_split = seen_resolved_clades[cur_parent][1]
        output[cur_parent] = (cur_split[1], cur_split[2])
        if len(cur_split[1]) > 1:
            working_list.append(cur_split[1])
        if len(cur_split[2]) > 1:
            working_list.append(cur_split[2])
    # reccursive function to write newick based on the output dict created above

    def recursive_nwk_split_dict(clade):
        nonlocal output
        nonlocal blockcount_map
        nonlocal branch_lengths_map
        if len(clade) == 1:
            return (f"{next(iter(clade.clade))}"
                    f"[&blockcount={np.median(blockcount_map[clade]) if clade.has_block else -1}]"
                    f":{np.mean(branch_lengths_map[clade])}")

        return (f"({recursive_nwk_split_dict(output[clade][0])},"
                f"{recursive_nwk_split_dict(output[clade][1])})"
                f"[&blockcount={np.median(blockcount_map[clade]) if clade.has_block else -1}]"
                f":{np.mean(branch_lengths_map[clade])}")

    return recursive_nwk_split_dict(root_clade)


def transmission_ccd_map_nexus(input_trees_file, output_tree_file, overwrite=False, burnin=0):
    """Takes input trees file and writes a tCCD1-MAP tree to the output file in Nexus format."""
    if os.path.exists(output_tree_file):
        if not overwrite:
            raise FileExistsError("File exists, set 'overwrite=True' to enable overwriting.")
        os.remove(output_tree_file)

    if not os.path.exists(input_trees_file):
        raise FileNotFoundError(f"Input trees {input_trees_file} not found.")

    if not 0 <= burnin < 1:
        raise ValueError("Burnin should be a number between 0 and 1, representing the proportion "
                         "of trees to delete at the beginning of the file.")

    trees = read_transmission_nexus_history(input_trees_file)
    trees = trees[int(0.1 * len(trees)):]  # deleting burnin from trees
    if len(trees) < 1:
        raise ValueError("Treeset is empty, reduce burning or check file.")

    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees)
    newick_map = get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)

    with (open(input_trees_file, 'r', encoding="UTF-8") as infile,
          open(output_tree_file, 'w+', encoding="UTF-8") as outfile):
        for line in infile:
            if line.strip().startswith("tree "):
                break  # Stop reading when the\ tree section starts
            outfile.write(line)

        outfile.write(f"tree tCCD_MAP = {newick_map};\nEnd;\n")
