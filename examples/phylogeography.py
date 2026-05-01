from pathlib import Path

from brokilon.ccd.domain.phylogeography import get_geo_map
from brokilon.core import read_nexus_trees


def read_phygeo_nexus():
    test_tree_file = (f"{Path(__file__).parent.absolute().parent}/"
                      f"examples/data/h3n2-bdmm.h3n2_2deme.trees")
    trees = read_nexus_trees(test_tree_file)
    # print(len(trees))
    # print(trees[123].write(features=["type"], format_root_node=True))

    geo_ccd_map, branch_lengths_map, clade_count_map = (
        get_geo_map(trees, geo_ann_str="type", ccd_type=1)
    )

    from brokilon.ccd.domain.phylogeography import sample_trees_from_geo_ccd
    random_ccd_sample = sample_trees_from_geo_ccd(10, geo_ccd_map, "type", clade_count_map)
    print(len(random_ccd_sample))
    for t in random_ccd_sample:
        # just for testing, this prints plain newick string now annotation...
        print(t.write())


if __name__ == '__main__':
    read_phygeo_nexus()
