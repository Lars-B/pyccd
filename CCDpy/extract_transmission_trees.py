import re
from CCDpy.label_transmission_history import label_transmission_tree

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
                # unknown_count = 0  # counting the unknown ancestors for separation


                # todo maybe make an unknown/transm_ancest class that we can easily check for instead of just a string...

                label_transmission_tree(new_tree)

                # rec_transmission_labelling(new_tree)

                trees.append(new_tree)
    return trees

