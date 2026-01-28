import os
from pathlib import Path

from brokilon.ccd import read_breath_nexus
from brokilon.ccd.domain.phylogeography import get_geo_map, get_geo_map_tree
from brokilon.core import read_nexus_trees


def phylgeo_example():
    posterior = read_nexus_trees(
        os.path.join(
            Path(__file__).parent.absolute(),
            "data/motivating_example.trees"
        )
    )

    print(len(posterior))
    # counts = {1: 4, 2: 2, 3: 2, 4: 1, 5: 1}
    # posterior = [trees[i] for i in range(len(trees)) for k in range(counts[i+1])]

    # todo think about ccd0 too... is there a difference?
    geo_ccd_map, branch_lengths_map, clade_count_map = (
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
    print(geo_ccd_map_tree.write(format=5, features=["type"],
                                 format_root_node=True))

    # todo get the regular CCD and map tree too
    # call beast for this to have the summary of states too...
    beast_ccd1_tree = read_nexus_trees(
        os.path.join(
            Path(__file__).parent.absolute(),
            "data/motivating_ccd1_beast.tree"
        )
    )[0]
    print("topo CCD1 map:")
    print(beast_ccd1_tree.write(format=5, features=["type"],
                                format_root_node=True))

    # todo can I reimplement the function to get all represented trees from a CCD?
    # todo implement the sample form the CCD function....
    from brokilon.ccd.domain.phylogeography import get_all_trees_represented
    all_trees_represented = get_all_trees_represented(geo_ccd_map)


def breath_example(basic=False):
    trees_path = (
        "data/basic_motivating_example_breath.trees" if basic else
        "data/motivating_example_breath.trees"
    )
    posterior = read_breath_nexus(
        os.path.join(
            Path(__file__).parent.absolute(),
            trees_path
        ),
        parse_taxon_map=False
    )

    print(len(posterior))

    from brokilon.ccd.domain.transmission import get_transmission_maps, \
        get_transmission_map_tree

    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(
        posterior,
        type_str="Ancestry")
    newick_map = get_transmission_map_tree(
        m1,
        m2,
        blockcount_map,
        branch_lengths_map,
        seed=1337
    )

    print("Transmission CCD1 map:")
    print(newick_map)
    print("-------")

    beast_ccd_tree_file = (
        "data/basic_motivating_breath_ccd1_beast.tree" if basic else
        "data/motivating_breath_ccd1_beast.tree"
    )
    beast_ccd1_tree = read_breath_nexus(
        os.path.join(
            Path(__file__).parent.absolute(),
            beast_ccd_tree_file
        ),
        False
    )[0]
    print("topo CCD1 map:")
    print(
        beast_ccd1_tree.write(
            format=5,
            features=["blockcount", "blockcount_median"],
            format_root_node=True
        )
    )


if __name__ == '__main__':
    phylgeo_example()
    # breath_example()
    # breath_example(basic=True)
