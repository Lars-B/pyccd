# todo sort this out as one thing,
#  how to use package data in here?
# def test_read_transmission_nexus():
#     """
#         WIP This will be moved to example folder shortly
#     """
#     test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
#     trees = pyccd.read_nexus.read_nexus_trees(test_tree_file, breath_trees=True,
#                                               label_transm_history=False)
#     m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees)
#     get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)
#     assert True, "Just testing if nothing breaks, apparently something did..."
#
#
# def test_transmisison_ccd_hiostory():
#     """
#         WIP This will be moved to example folder shortly
#     """
#     test_tree_file = f"{Path(__file__).parent.absolute()}/data/Filter-roetzer40.trees"
#     trees = pyccd.read_nexus.read_nexus_trees(test_tree_file, breath_trees=True,
#                                               label_transm_history=True)
#     m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees, type_str="Ancestry")
#     tree = get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)
#     assert tree is not None, "Failed this example..."
