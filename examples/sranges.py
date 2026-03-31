from pathlib import Path


test_tree_file = (f"{Path(__file__).parent.absolute().parent}"
                  f"/examples/data/sr_example.trees")

test_mcc_tree_file = (f"{Path(__file__).parent.absolute().parent}"
                  f"/examples/data/sr_mcc.tree")

from brokilon.core import read_nexus_trees

trees, map = read_nexus_trees(test_tree_file, parse_taxon_map=True)

from brokilon.ccd.domain.sranges import sranges

clade_counts, clade_split_counts = sranges.get_sranges_map(trees, map)

taxon_map = map
reverse_taxon_map = {value: key for key, value in taxon_map.items()}

map_tree = sranges.get_sranges_map_tree(clade_counts, clade_split_counts, taxon_map, reverse_taxon_map)


# for c in clade_counts.keys():
#     print(c)
#
# for c in clade_split_counts.keys():
#     print(c)
# print(clade_counts)
# print(clade_split_counts)
