from pathlib import Path

from brokilon import read_nexus_trees


def read_phygeo_nexus():
    test_tree_file = (f"{Path(__file__).parent.absolute().parent}/"
                      f"examples/data/h3n2-bdmm.h3n2_2deme.trees")
    trees = read_nexus_trees(test_tree_file, breath_trees=False)
    print(len(trees))
    print(trees[123].write(features=["type"], format_root_node=True))


if __name__ == '__main__':
    read_phygeo_nexus()
