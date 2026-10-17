from pathlib import Path


# test_tree_file = (f"{Path(__file__).parent.absolute().parent}"
#                   f"/examples/data/sr_example.trees")
# test_tree_file = (f"{Path(__file__).parent.absolute().parent}"
#                   f"/examples/data/sr_small.trees")

test_tree_file = (f"{Path(__file__).parent.absolute().parent}"
                  f"/examples/data/sranges_toy_posterior.trees")

# test_tree_file = (f"{Path(__file__).parent.absolute().parent}"
#                   f"/examples/data/rep_3_srfbd_first_ucln.trees")

# test_mcc_tree_file = (f"{Path(__file__).parent.absolute().parent}"
#                   f"/examples/data/sr_mcc.tree")

from brokilon.core import read_nexus_trees

trees, map = read_nexus_trees(test_tree_file, parse_taxon_map=True)

# subset for debugging
# single range birth
# trees = trees[11:12]
# stacked ranges case
# trees = trees[26:27]

# trees = trees[6:7]

print(f"Looking at {len(trees)} many trees.")
from brokilon.ccd.domain.sranges import sranges

clade_counts, clade_split_counts, sranges_set, sampled_ancestors = sranges.get_sranges_map(trees, map)

taxon_map = map

reverse_taxon_map = {value: key for key, value in taxon_map.items()}

# making a dot graph of the CCD

from graph_generation import sranges_map_to_networkx

# G = sranges_map_to_networkx(clade_split_counts)

from brokilon.ccd.domain.sranges.sranges_ccd_graph_visualizer import sranges_map_to_cytoscape_html

sranges_map_to_cytoscape_html(
    clade_split_counts, reverse_taxon_map, clade_counts,
    use_reverse_taxon_map=False
)

print(f"The following are ranges:\n{sranges_set}")

sranges.DEBUG_TREE = False

seen_resolved, map_tree = sranges.get_sranges_map_tree(
        clade_counts,
        clade_split_counts,
        sranges_set,
        sampled_ancestors
)

nwk_map = map_tree.write(
    format=1,
    format_root_node=True,
    features=["orientation", "ancestral_range"],
    icytree_order=True
)
print(nwk_map)

# for i in range(len(trees)):
#     clade_counts, clade_split_counts = sranges.get_sranges_map([trees[i]], map)
#
#     taxon_map = map
#     reverse_taxon_map = {value: key for key, value in taxon_map.items()}
#
#     testing = [[k for k in clade_counts if k.clade == c.clade] for c in clade_counts]
#
#     testing_concrete = [l for l in testing if len(l) > 1]
#
#     for l in testing:
#         if len(l) != 1:
#             print("This is a problem we need to fix")
#
#     map_tree = sranges.get_sranges_map_tree(
#         clade_counts,
#         clade_split_counts,
#         taxon_map,
#         reverse_taxon_map
#     )
#     print(f"{i}. Finished...")

# for c in clade_counts.keys():
#     print(c)
#
# for c in clade_split_counts.keys():
#     print(c)
# print(clade_counts)
# print(clade_split_counts)
