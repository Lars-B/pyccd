import os
from pathlib import Path

from brokilon.ccd.domain.phylogeography import get_geo_map, get_geo_map_tree
from brokilon.core import read_nexus_trees


if __name__ == '__main__':
    trees = read_nexus_trees(
        os.path.join(
            Path(__file__).parent.absolute(),
            "data/motivating_example.trees"
        )
    )

    print(len(trees))
    counts = {1: 4, 2: 2, 3: 2, 4: 1, 5: 1}
    posterior = [trees[i] for i in range(len(trees)) for k in range(counts[i+1])]

    # todo think about ccd0 too... is there a difference?
    geo_ccd_map, branch_lenghts_map, clade_count_map = (
        get_geo_map(posterior, geo_ann_str="type", ccd_type=1)
    )
    geo_ccd_map_tree = (
        get_geo_map_tree(
            geo_ccd_map,
            geo_ann_str="type",
            clade_count_map=clade_count_map
        )
    )
    print("Geo CCD1 map:")
    print(geo_ccd_map_tree.write(format=5, features=["type"], format_root_node=True))

    # todo get the regular CCD and map tree too
    # call beast for this to have the summary of states too...
    beast_ccd1_tree = read_nexus_trees(
        os.path.join(
            Path(__file__).parent.absolute(),
            "data/motivating_ccd1_beast.tree"
        )
    )[0]
    print("topo CCD1 map:")
    print(beast_ccd1_tree.write(format=5, features=["type"], format_root_node=True))
    # TODO Note that I currently have no idea how beast gets the root annotation? why is it 0+1
    #  and why are the probabilities of 0 and 1 0.4?

    # todo can I reimplement the function to get all represented trees from a CCD?
    # todo implement the sample form the CCD function....

