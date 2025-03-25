import re
from logging import CRITICAL

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
                # tree_string = f'{re.split(re_tree, line)[1][:re.split(re_tree, line)[1].rfind(")") + 1]};'
                tree_string = re.split(re_tree, line)[1]
                pattern = r"\d*\[[^\]]*\]"

                # matches = re.findall(pattern, tree_string)

                counter = 0  # Initialize a counter

                def replace_match(match):
                    nonlocal counter
                    split_str = match.group(0).split("[&")
                    if len(split_str) == 2:
                        meta_list = split_str[1][:-1].split(
                            ",")  # deleting the ] from string and splitting all the meta data
                        matching_element = next((s for s in meta_list if "blockcount" in s), None)
                        block_count = int(
                            float(matching_element.split("=")[-1]))  # variable hard code with string in line above
                        if not split_str[0]:
                            counter += 1
                        # todo add another / and the ancestral infector
                        return f"%{split_str[0] if split_str[0] else f'internal{counter}'}/{block_count}%"  # making new string be pair taxa;block_count
                    else:
                        raise NotImplementedError("Needs to be added or debugged.")

                # Replace all matches
                new_tree_string = re.sub(pattern, replace_match, tree_string)
                new_tree = ete3.Tree(new_tree_string, format=1)

                # taxa_labels = [l.name.split("/")[0][-1:] for l in new_tree]  # get a list of all taxa labels
                # for node in new_tree.traverse("levelorder"):
                #     if node.is_root():
                #         # Root as special case
                #         transm_map[node.name] = "unknown"
                #     else:
                #         # non root vertex'
                #         cur_name, cur_blockcount = node.name.replace("%", "").split("/")
                #         cur_parent_name = node.up.name
                #         if cur_parent_name not in transm_map:
                #             raise AttributeError("Traversal through tree failed, parent should be in the transmission map at this point!")
                #         if cur_blockcount == -1:
                #             # the edge towards the current node had not transmission event
                # # todo traverse the new tree to label the transmission history into the nodes
                unknown_count = 0  # counting the unknown ancestors for separation
                transmission_ancestor_key = "transm_ancest"

                # todo maybe make an unknown/transm_ancest class that we can easily check for instead of just a string...

                def rec_transmission_labelling(node):
                    nonlocal unknown_count
                    nonlocal transmission_ancestor_key
                    if "blockcount" in node.features and transmission_ancestor_key in node.features:
                        # We have visited this node and labeled it, nothing to do with this node
                        # todo this might cause premature termination??? Second part of if...
                        return
                    if "blockcount" not in node.features:
                        # chekc that we haven't updated yet because the code might visit nodes multiple times
                        # updating the node name and blockcount information
                        cur_name, cur_blockcount = node.name.replace("%", "").split("/")
                        node.add_feature("blockcount", int(cur_blockcount))
                        node.name = cur_name

                    cur_unknonw = "Currently Unknown"

                    # Todo these are the possible cases for a clade
                    #  We have a block coming into the node therefore the ancestor is UNKNOWN
                    #  We have a single event coming in to the node, the ancestor is either UNKNOWN or any SAMPLE
                    #  We have no event coming into the node, the ancestor is the same as the one of the node above
                    #  try to only ever label downward

                    if node.is_leaf():
                        if node.blockcount == -1:
                            node.add_feature(transmission_ancestor_key, node.name)
                        elif node.blockcount == 0:
                            node.add_feature(transmission_ancestor_key, f"Unknown-{unknown_count}")
                            unknown_count += 1
                        elif node.blockcount > 0:  # todo this can be differentiated into block and no block, in case of block we just set it to unknown!
                            # get the other sibling of current node
                            if node.up.children[0] != node:
                                other_sibling_node = node.up.children[0]
                            else:
                                other_sibling_node = node.up.children[1]
                            # if the subtree rooted at the other sibling is not labeled, call the recursion
                            node.add_feature(transmission_ancestor_key, cur_unknonw)  # to stop infinite recursion below
                            if not transmission_ancestor_key in other_sibling_node.features:
                                rec_transmission_labelling(other_sibling_node)
                            # now we know that the other_sibling_node has a label transm_hist
                            if (other_sibling_node.blockcount) >= 0:
                                # This is the case where both the current leaf edge and the edge to the other sibling have an event
                                # Therefore, we need to look at the current parent of both nodes
                                if node.up.blockcount >= 0:
                                    # event before this node.up and afterwards towards both children: Truly unknown sample
                                    # The history for all three nodes in question is an unknown host
                                    node.up.transm_ancest = f"Unknown-{unknown_count}"
                                    node.transm_ancest = f"Unknown-{unknown_count}"
                                    other_sibling_node.transm_ancest = f"Unknown-{unknown_count}"
                                    unknown_count += 1
                                elif node.up.blockcount == -1:
                                    # todo this is weird and currently not necessary...
                                    # todo add a tree to the test case where this happens
                                    # todo check if this is even allowed...
                                    print("bla")
                                else:
                                    raise ValueError("The blockcount should be -1, 0, or some positive integer!")
                            elif int(other_sibling_node.blockcount) == -1:
                                # No mutation of the other sibling, therefore the ancestor of current parent and current node is the same as that one
                                node.transm_ancest = other_sibling_node.transm_ancest
                        else:
                            raise ValueError("Blockcount should not be less than -1.")
                    else:
                        child0 = node.children[0]
                        child1 = node.children[1]
                        rec_transmission_labelling(child0)
                        rec_transmission_labelling(child1)
                        if transmission_ancestor_key in node.features:
                            # todo make sure this actually doesn't happen once finished
                            raise ValueError("I don't think this case should happen, maybe...?")
                        else:
                            if child0.blockcount >= 0 and child1.blockcount >= 0:
                                if node.is_root() or node.blockcount >= 0:
                                    # Either root or in a state where there is a transmission event above and below towards both children
                                    # i.e. an unknown sample
                                    node.transm_ancest = f"Unknown-{unknown_count}"
                                    child0.transm_ancest = f"Unknown-{unknown_count}"
                                    child1.transm_ancest = f"Unknown-{unknown_count}"
                                    unknown_count += 1
                                else:
                                    # Case transmission event towards both children, but not towards current node
                                    print("This is a stupid case")
                                    # todo do same as for leafs, move one step up, label sibling, check again...
                            if child0.blockcount == child1.blockcount == -1:
                                # Internal node where nothing happens towards both children
                                # In this case need to look further up the tree for answers
                                print("This is also a stupid case...")

                            # Here we know either of the two children has blockcount == -1
                            child0, child1 = (child0, child1) if child0.blockcount == -1 else (child1, child0)
                            # Child0 has nothing happen towards it, but child 1 does
                            if not child1.blockcount >= 0:
                                raise ValueError(
                                    "This should definitely not happen! Substantial problem within tree and/or code.")
                            if node.blockcount == 0:
                                # Case:
                                # child0 has no event towards it
                                # child1 has an event towards
                                # The current node has a single event towards it
                                # todo this case is currently happening...
                                print("one event, therefore we might know the history")
                            elif node.blockcount > 0:
                                # Case:
                                # child0 has no event towards it
                                # child1 has an event towards
                                # The current node has a block of events towards it
                                node.transm_ancest = f"Unknown-{unknown_count}"
                                child0.transm_ancest = f"Unknown-{unknown_count}"
                                child1.transm_ancest = f"Unknown-{unknown_count}"
                                unknown_count += 1
                            else:
                                # Case:
                                # child0 has no event towards it
                                # child1 has an event towards
                                # The current node has no event towards it
                                if child0.transm_ancest == cur_unknonw:
                                    print("Another case where we need to look upwards again.")
                                else:
                                    node.transm_ancest = child0.transm_ancest
                                    if child1.blockcount == 0:
                                        child1.transm_ancest = node.transm_ancest

                rec_transmission_labelling(new_tree)

                trees.append(new_tree)
                raise ValueError("One tree done for now.")
    return trees
