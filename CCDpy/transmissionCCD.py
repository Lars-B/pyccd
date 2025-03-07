import re
from collections import defaultdict
from dataclasses import dataclass

import ete3


def read_transmission_nexus(file: str, _ignore_ctrees=False) -> list:
    # re_tree returns nwk string without the root height and no ; in the end
    re_tree = re.compile("\t?tree .*=? (.*$)", flags=re.I | re.MULTILINE)
    # Used to delete the ; and a potential branchlength of the root
    # name_dict = get_mapping_dict(file)  # Save tree label names in dict

    # todo this will be changed once we are interested in the branch annotation information
    brackets = re.compile(r'\[[^\]]*\]')  # Used to delete info in []

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


def get_transmission_clades(tree, blockcountmap):
    treestr = tree.write(format=8)  # todo this ignores the internal node names, which is wrong
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
                    blockcountmap[curclade] += blockcount
                else:
                    blockcountmap[curclade] = blockcount
            clades.add(curclade)
            del opend[-1]
    if opend:
        raise Exception("Invalid tree string given! (to many '(')")
    return clades


def get_transmission_maps(trees):

    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled
    uniques = {}

    seen = {}
    # todo m1 needs to have the root and the single leaf clades
    #  need to get a list of weights for the unique trees
    #  this needs to replace trees in the for loop and the weight will be added instead of 1 each time!

    for ix, t in enumerate(trees):
        if not frozenset(sorted(get_transmission_clades(t))) in seen:
            seen[frozenset(sorted(get_transmission_clades(t)))] = ix
            uniques[ix] = []
        else:
            uniques[seen[frozenset(sorted(get_transmission_clades(t)))]].append(ix)

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