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

trees, taxon_map = read_nexus_trees(test_tree_file, parse_taxon_map=True)

# subset for debugging
# single range birth
# trees = trees[11:12]
# stacked ranges case
# trees = trees[26:27]
# trees = trees[6:7]

print(f"Looking at {len(trees)} trees.")
from brokilon.ccd.domain.sranges import sranges

# todo this can only be called once, the prelabel thing needs to be done better
#  i need to refactor that code properly...
# clade_counts, clade_split_counts, sranges_set, sampled_ancestors = sranges.get_sranges_map(trees, taxon_map)

reverse_taxon_map = {value: key for key, value in taxon_map.items()}

from brokilon.metrics import extended_robinson_foulds

# dist = extended_robinson_foulds(trees[0], trees[1], kind="sranges", taxon_map=taxon_map)

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

n = len(trees)
dist_matrix = np.zeros((n, n), dtype=int)

for i, j in combinations(range(n), 2):
    dist = extended_robinson_foulds(
        trees[i], trees[j], kind="sranges", taxon_map=taxon_map
    )
    dist_matrix[i, j] = dist
    dist_matrix[j, i] = dist

plt.figure(figsize=(max(8, n * 0.6), max(6, n * 0.5)))

sns.heatmap(
    dist_matrix,
    annot=True,
    fmt="d",
    cmap="viridis",
    xticklabels=range(n),
    yticklabels=range(n),
    cbar=False,
)

plt.xlabel("Tree")
plt.ylabel("Tree")
plt.title("Pairwise Extended Robinson–Foulds Distances")
plt.tight_layout()
plt.savefig("heatmap_eRF_splits.pdf")


# making a dot graph of the CCD

from graph_generation import sranges_map_to_networkx

# G = sranges_map_to_networkx(clade_split_counts)

# from brokilon.ccd.domain.sranges.sranges_ccd_graph_visualizer import sranges_map_to_cytoscape_html

# sranges_map_to_cytoscape_html(
#     clade_split_counts, reverse_taxon_map, clade_counts,
#     use_reverse_taxon_map=False
# )

# print(f"The following are ranges:\n{sranges_set}")

# sranges.DEBUG_TREE = False

# seen_resolved, map_tree = sranges.get_sranges_map_tree(
#         clade_counts,
#         clade_split_counts,
#         sranges_set,
#         sampled_ancestors
# )
#
# nwk_map = map_tree.write(
#     format=1,
#     format_root_node=True,
#     features=["orientation", "ancestral_range"],
#     icytree_order=True
# )
# print(nwk_map)

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
