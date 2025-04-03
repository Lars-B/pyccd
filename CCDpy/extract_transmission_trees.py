"""
This module reads transmission nexus files as returned by the BREATH BEAST2 package.
It will read the trees into ete3.Tree objects and label them with their transmission ancestry.

Dependencies:
- ete3.Tree
"""
import re
import ete3

from CCDpy.label_transmission_history import label_transmission_tree



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
    # Used to delete the ; and a potential branch length of the root
    # name_dict = get_mapping_dict(file)  # Save tree label names in dict

    trees = []
    with open(file, 'r', encoding="utf-8") as f:
        for line in f:
            if re_tree.match(line):
                # tree_string = f'{re.split(re_tree, line)[1][:re.split(re_tree, line)[1]
                #                                              .rfind(")") + 1]};'
                tree_string = re.split(re_tree, line)[1]
                pattern = r"\d*\[[^\]]*\]"

                counter = 0  # Initialize a counter

                def replace_match(match):
                    nonlocal counter
                    split_str = match.group(0).split("[&")
                    if len(split_str) == 2:
                        meta_list = split_str[1][:-1].split(
                            ",")  # deleting the ] from string and splitting all the meta data
                        matching_element = next((s for s in meta_list if "blockcount" in s), None)
                        # variable hard code with string in line above
                        block_count = int(
                            float(matching_element.split("=")[-1]))
                        if not split_str[0]:
                            counter += 1
                        # making new string be pair taxa;block_count
                        node_name = split_str[0] if split_str[0] else f'internal{counter}'
                        return f"%{node_name}/{block_count}%"
                    raise NotImplementedError("Needs to be added or debugged.")

                # Replace all matches
                new_tree_string = re.sub(pattern, replace_match, tree_string)
                new_tree = ete3.Tree(new_tree_string, format=1)


                # adjusting the tree to contain the blockcount label and correct node names
                for node in new_tree.traverse("levelorder"):
                    # this should technically never be the case...
                    if not hasattr(node, "blockcount"):
                        # Assert node.name format
                        assert node.name and "/" in node.name, (f"Invalid node name format: "
                                                                f"'{node.name}' "
                                                                f"(expected 'name/blockcount')")

                        # Extract node name and blockcount safely
                        # Ensure only one split
                        cur_name, cur_blockcount = node.name.replace("%", "").split("/", 1)
                        # Strip spaces & convert safely
                        node.add_feature("blockcount", int(cur_blockcount.strip()))
                        # Remove unnecessary spaces
                        node.name = cur_name.strip()

                label_transmission_tree(new_tree)

                # rec_transmission_labelling(new_tree)

                trees.append(new_tree)
    return trees
