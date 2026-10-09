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

# trees = trees[:10]

from brokilon.ccd.domain.sranges import sranges

clade_counts, clade_split_counts, sranges_set, sampled_ancestors = sranges.get_sranges_map(trees, map)

taxon_map = map

reverse_taxon_map = {value: key for key, value in taxon_map.items()}

# making a dot graph of the CCD

from graph_generation import sranges_map_to_networkx

# G = sranges_map_to_networkx(clade_split_counts)

from ccd_vis import sranges_map_to_cytoscape_html

# todo next steps to get this going:
#  - improve visualization with probabilities on edges
#  - need to make sure that we are getting the
#    right maps of clades and splits....

sranges_map_to_cytoscape_html(clade_split_counts, reverse_taxon_map, clade_counts)

seen_resolved, map_tree = sranges.get_sranges_map_tree(
        clade_counts,
        clade_split_counts,
        sranges_set,
        sampled_ancestors,
        taxon_map,
        reverse_taxon_map
)

print(f"length of seen resolved: {len(seen_resolved.keys())}")
for clade in seen_resolved:
    cur_c = str({reverse_taxon_map[c] for c in clade.clade})
    cur_range = reverse_taxon_map[clade.ancestral_range] if clade.ancestral_range is not None else "NR"
    out = f"{cur_c}_{cur_range}"
    print(out)

nwk_map = map_tree.write(
    format=1,
    format_root_node=True,
    features=["orientation", "ancestral_range"]
)

# parent_node.up.up.write(format=1, features=["orientation", "ancestral_range"])

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
