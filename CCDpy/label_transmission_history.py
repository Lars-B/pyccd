def label_transmission_tree(etree):
    # separately label the nodes with blockcount, this should be included in the thing below later on...
    for node in etree.traverse("levelorder"):
        if "blockcount" not in node.features:
            # chekc that we haven't updated yet because the code might visit nodes multiple times
            # updating the node name and blockcount information
            cur_name, cur_blockcount = node.name.replace("%", "").split("/")
            node.add_feature("blockcount", int(cur_blockcount))
            node.name = cur_name

    unknown_count = 0

    unlabeled_nodes_list = []
    top_infected_nodes_list = []

    # Label all nodes that have a block before them, these have unknonw ancestry
    for node in etree.traverse("levelorder"):
        # Labeled nodes don't have to be considered
        if hasattr(node, "transm_ancest"):
            continue
        # There can still be leaves here that have an unknown transmission ancestor...
        if node.blockcount > 0:
            # There is a block so we can immediately lable this transmission ancestry
            node.transm_ancest = f"Unknown-{unknown_count}"
            unknown_count += 1
        else:
            unlabeled_nodes_list.append(node)

    # First label all leafs and reachable nodes from leaves
    for leaf in etree:
        if leaf.blockcount == -1:
            # No transmission event towards the current node, recursively label all the nodes reachable from here
            leaf.transm_ancest = leaf.name
            if leaf in unlabeled_nodes_list:
                unlabeled_nodes_list.remove(leaf)
            working_node_list = [leaf.up]

            while working_node_list:
                working_node = working_node_list.pop()
                if working_node.blockcount == -1:  # and not working_node.is_root():
                    # THere is no event towards current node
                    working_node.transm_ancest = leaf.name  # label the current node
                    if working_node in unlabeled_nodes_list:
                        unlabeled_nodes_list.remove(working_node)
                    # if the current node is not the root we can look further up
                    if working_node.up:
                        if not hasattr(working_node.up, "transm_ancest"):
                            # if the node upwards has not been labeled we can check if we need to add it
                            if working_node.up.blockcount == -1:
                                # the node up has not been labeled and extends the path
                                if not working_node.up in working_node_list:
                                    working_node_list.append(working_node.up)

                if len(working_node.children) > 0:
                    assert len(working_node.children) == 2, "Non binary tree found, not supported!"
                    for c in working_node.children:
                        # adding the children that can be reached from current node
                        if not hasattr(c, "transm_ancest"):
                            if c.blockcount == -1:
                                # adding the children that are reachable but have not been labeled yet
                                if not c in working_node_list:
                                    working_node_list.append(c)
                            elif c.blockcount == 0:
                                # child that is infected by current label
                                c.transm_ancest = leaf.name
                                top_infected_nodes_list.append(c)
                                if c in unlabeled_nodes_list:
                                    # Removing the child c if it was previously added to unlabeled list
                                    unlabeled_nodes_list.remove(c)
        elif leaf.blockcount == 0 and leaf not in top_infected_nodes_list and leaf not in unlabeled_nodes_list:
            unlabeled_nodes_list.append(leaf)

    while top_infected_nodes_list:
        cur_node = top_infected_nodes_list.pop()
        assert hasattr(cur_node, "transm_ancest"), "This node has to be labeled here!"
        if cur_node in unlabeled_nodes_list:
            print("This shouldn't happen here...!")
            unlabeled_nodes_list.remove(cur_node)
        if len(cur_node.children) > 0:
            assert len(cur_node.children) == 2, "Non binary tree found, not supported!"
            for c in cur_node.children:
                if not hasattr(c, "transm_ancest"):
                    if c.blockcount in (0, -1):
                        # child of top infected node
                        # if no event or one event on the edge we can label it with the ancestor
                        c.transm_ancest = cur_node.transm_ancest
                        if c.blockcount == -1:
                            # only need to keep it if there was no event, i.e. path continues
                            top_infected_nodes_list.append(c)
                            if c in unlabeled_nodes_list:
                                # Removing the child c if it was previously added to unlabeled list
                                unlabeled_nodes_list.remove(c)


    if len(unlabeled_nodes_list) == 0:
        return

    # todo might want to use this for etree: t.search_nodes(attr=value)

    # Sorting the list of unlabeled nodes by their level to do a for loop that is doing level order within the whole tree
    unlabeled_nodes_list.sort(key=lambda obj: obj.get_distance(etree.get_tree_root(), topology_only=True))
    for node in unlabeled_nodes_list:
        assert not hasattr(node, "transm_ancest"), "This should never happen! We visited a node twice for no need..."
        assert node.blockcount < 1, "These nodes should be labeled at this point."
        # all the paths from leaves have been labeled, at this point there is no path to any leaf from here
        # since we have ordered the nodes in the list accordign to their 'level' in the tree the ancestry pulls down
        if node.is_root():
            node.transm_ancest = f"Unknown-{unknown_count}"
            unknown_count += 1
        elif node.up.blockcount == -1:
            node.transm_ancest = node.up.transm_ancest
        else:
            node_sibling = next(x for x in node.up.children if x != node)
            if node_sibling.blockcount == -1:
                node.transm_ancest = node_sibling.transm_ancest
            else:
                node.transm_ancest = f"Unknown-{unknown_count}"
                unknown_count += 1

    # todo add a counter of how often we visit each node, try to minimize that to be 1 at best
    # todo write the BREATH testing tree file so that we have all the cases visited at least once...
    cur_tree_nwk = etree.write(features=["transm_ancest"], format_root_node=True, format=2)
    return