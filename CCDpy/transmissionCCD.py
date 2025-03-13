import os.path
import re
from collections import defaultdict
from dataclasses import dataclass
from idlelib.pyparse import trans
from multiprocessing.managers import Value

import ete3


def read_transmission_nexus(file: str) -> list:
    """
    Function to read a nexus file that contains transmission trees.

    :param file: Input file
    :type file: str
    :return: list of transmission trees
    :rtype: list
    """
    # re_tree returns nwk string without the root height and no ; in the end
    re_tree = re.compile("\t?tree .*=? (.*$)", flags=re.I | re.MULTILINE)
    # Used to delete the ; and a potential branchlength of the root
    # name_dict = get_mapping_dict(file)  # Save tree label names in dict

    trees = []
    with open(file, 'r') as f:
        for line in f:
            if re_tree.match(line):
                tree_string = f'{re.split(re_tree, line)[1][:re.split(re_tree, line)[1].rfind(")") + 1]};'
                # todo create a fully labelled tree, i.e. each node is labelled with a tuple
                #  the tuple contains (taxa id, block count)
                #  this should be able to contain more than just a specific information?

                pattern = r"\d*\[[^\]]*\]"

                # matches = re.findall(pattern, tree_string)

                counter = 0  # Initialize a counter
                def replace_match(match):
                    nonlocal counter
                    counter += 1
                    split_str = match.group(0).split("[&")
                    if len(split_str) == 2:
                        # its a leaf edge so there is already a taxon name
                        meta_list = split_str[1][:-1].split(",") # deleting the ] from string and splitting all the meta data
                        matching_element = next((s for s in meta_list if "blockcount" in s), None)
                        block_count = int(
                            float(matching_element.split("=")[-1]))  # variable hard code with string in line above
                        return f"%{split_str[0]}/{block_count}%"  # making new string be pair taxa;block_count
                    elif len(split_str) == 1:
                        # internal edge so no node label
                        meta_list = split_str[0][:-1].split(",")  # deleting the ] from string and splitting all the meta data
                        matching_element = next((s for s in meta_list if "blockcount" in s), None)
                        block_count = int(
                            float(matching_element.split("=")[-1]))  # variable hard code with string in line above
                        # todo this naming doesn't work properly
                        return f"%internal{counter}/{block_count}%"  # making new string be pair taxa;block_count
                    else:
                        raise NotImplementedError("Needs to be added or debugged.")
                # Replace all matches
                new_tree_string = re.sub(pattern, replace_match, tree_string)
                trees.append(ete3.Tree(new_tree_string, format=1))
                # trees.append(ete3.Tree(re.sub(brackets, "", tree_string), format=0))
    return trees

@dataclass(frozen=True)
class transmission_clade:
    '''Class to keep clades with additional transmission block information'''
    clade: frozenset
    has_block: bool
    # blockcount: int  # not sure if I want this in here or not...
    def __len__(self):
        return len(self.clade)

    def __lt__(self, other):
        if isinstance(other, transmission_clade):
            return len(self.clade) < len(other.clade)
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, transmission_clade):
            return len(self.clade) > len(other.clade)
        return NotImplemented


def get_transmission_clades(tree, blockcountmap):
    treestr = tree.write(format=8)
    clades = set()
    if not (treestr[0] == '(' and treestr[-2] == ')' and treestr[-1] == ';'):
        # todo the above second statement is not true since there might be root info that we can ignore
        raise Exception("Invalid tree string given! (no ';' at the end)")

    # Adding leafs separately
    leafnamepattern = r"%(-?\d+/-?\d+)%"
    # the above regex assumes that the taxa are integers and internal nodes do not have an interger in the first part of the block.
    matches = re.findall(leafnamepattern, treestr)
    for m in matches:
        taxa, blockcount = map(int, m.split("/"))
        if blockcount == -1:
            leafclade = transmission_clade(frozenset({taxa}), False)
        else:
            leafclade = transmission_clade(frozenset({taxa}), True)
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
                raise Exception("Invalid tree string given! (to many ')')")
            # current_string = treestr[opend[-1]:i]

            match_interenalnodename = re.search(internalnodenamepattern, treestr[i+1:])
            if not match_interenalnodename:
                raise ValueError("Something went wrong looking for the next internal nodes name")
            blockcount = int(match_interenalnodename.group(0).replace("%", "").split("/")[1])
            current_clade_set = frozenset(re.sub(re_brackets_internals, "", treestr[opend[-1]:i]).split(","))
            if blockcount == -1:
                curclade = transmission_clade(current_clade_set, False)
            else:
                curclade = transmission_clade(current_clade_set, True)
                if curclade in blockcountmap:
                    blockcountmap[curclade].append(blockcount)
                else:
                    blockcountmap[curclade] = [blockcount]
            clades.add(curclade)
            del opend[-1]
    if opend:
        raise Exception("Invalid tree string given! (to many '(')")
    return clades


def get_transmission_maps(trees):

    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled
    blockcount_map = {}

    for ix, t in enumerate(trees):
        for node in t.traverse("levelorder"):
            # todo this if should be > 1 since we also care about leafs now.
            if len(node) > 1:
                c = node.children
                c0_leafs = set()
                for leaf in c[0]:
                    c0_leafs.add(int(leaf.name.replace("%", "").split("/")[0]))
                c1_leafs = set()
                for leaf in c[1]:
                    c1_leafs.add(int(leaf.name.replace("%", "").split("/")[0]))

                parent_clade_set = frozenset(sorted(c0_leafs.union(c1_leafs)))
                if node.name:
                    blockcount = int(node.name.replace("%", "").split("/")[1])
                    if blockcount == -1:
                        parent_clade = transmission_clade(parent_clade_set, False)
                    else:
                        parent_clade = transmission_clade(parent_clade_set, True)
                        if parent_clade in blockcount_map:
                            blockcount_map[parent_clade].append(blockcount)
                        else:
                            blockcount_map[parent_clade] = [blockcount]
                else:
                    # this is the root, apprarently we don't have the name atm
                    parent_clade = transmission_clade(parent_clade_set, False)

                m1[parent_clade] += 1

                blockcount_c0 = int(c[0].name.replace("%", "").split("/")[1])
                blockcount_c1 = int(c[1].name.replace("%", "").split("/")[1])
                child0_clade = transmission_clade(frozenset(c0_leafs), blockcount_c0 != -1)
                child1_clade = transmission_clade(frozenset(c1_leafs), blockcount_c1 != -1)

                # in m2 the split clade with the lower int for taxa is entered first
                if min(c0_leafs) < min(c1_leafs):
                    m2[(parent_clade, child0_clade, child1_clade)] += 1
                else:
                    m2[(parent_clade, child1_clade, child0_clade)] += 1
            elif len(node) == 1:
                # leaf node for which we need to add the blockcount to the blockcount_map
                leaf_label, blockcount = map(int, node.name.replace("%", "").split("/"))
                if blockcount != -1:
                    leaf_clade = transmission_clade(frozenset({leaf_label}), True)
                    if leaf_clade in blockcount_map:
                        blockcount_map[leaf_clade].append(blockcount)
                    else:
                        blockcount_map[leaf_clade] = [blockcount]
    return m1, m2, blockcount_map


def get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map):
    seen_resolved_clades = {}
    all_clades = sorted(list(m1.keys()), key=len)  # sorted list of all clades, small (cherries) to big

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
                # todo the following is not randomly resolving tie breaks, just keeps the first found
                if seen_resolved_clades[current_split[0]][0] <= split_prob:
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
    import numpy as np
    # reccursive function to write newick based on the output dict created above
    def recursive_nwk_split_dict(clade):
        nonlocal output
        nonlocal blockcount_map
        if len(clade) == 1:
            return f"{next(iter(clade.clade))}[&blockcount={np.median(blockcount_map[clade]) if clade.has_block else -1}]:1.0"
        else:
            return (f"({recursive_nwk_split_dict(output[clade][0])},"
                    f"{recursive_nwk_split_dict(output[clade][1])})"
                    f"[&blockcount={np.median(blockcount_map[clade]) if clade.has_block else -1}]:1.0")

    return recursive_nwk_split_dict(root_clade)


def transmissionCCD_MAP_nexus(input_trees_file, output_tree_file, overwrite=False, burnin=0):
    """Takes input trees file and writes a tCCD1-MAP tree to the output file in Nexus format."""
    if os.path.exists(output_tree_file):
        if not overwrite:
            raise FileExistsError("File exists, set 'overwrite=True' to enable overwriting.")
        os.remove(output_tree_file)

    if not os.path.exists(input_trees_file):
        raise FileNotFoundError(f"Input trees {input_trees_file} not found.")

    if not (0 <= burnin < 1):
        raise ValueError("Burnin should be a number between 0 and 1 representing the proportion of trees to delete at the beginning of the file.")

    trees = read_transmission_nexus(input_trees_file)
    trees = trees[int(0.1 * len(trees)):]  # deleting burnin from trees
    if len(trees) < 1:
        raise ValueError("Treeset is empty, reduce burning or check file.")

    m1, m2, blockcount_map = get_transmission_maps(trees)
    newick_MAP = get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map)

    #def copy_nexus_file(input_file, output_file):
    with open(input_trees_file, 'r') as infile, open(output_tree_file, 'w+') as outfile:
        for line in infile:
            if line.strip().startswith("tree "):
                break  # Stop reading when the\ tree section starts
            outfile.write(line)

        outfile.write(f"tree tCCD_MAP = {newick_MAP};\nEnd;\n")
