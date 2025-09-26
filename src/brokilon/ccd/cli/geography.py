import click
from pathlib import Path
from brokilon.core.read_nexus import read_nexus_trees
from brokilon.ccd.domain.phylogeography import get_geo_map, get_geo_map_tree
from brokilon.ccd.cli.common_options import common_options


@click.command()
@common_options
def main(trees_file, ccd_type, burn_in, output_file):
    """
    Process a Nexus tree file and generate a geo-mapped tree.
    """
    trees_file = Path(trees_file).absolute()

    # Read trees and taxon map
    trees, taxon_map = read_nexus_trees(trees_file, parse_taxon_map=True)

    # Skip a fraction of trees if needed
    start_idx = int(len(trees) * burn_in)
    trees = trees[start_idx:]

    # Compute geo map and branch lengths
    geo_ccd_map, branch_length_map = get_geo_map(trees, geo_ann_str="type", ccd_type=ccd_type)

    # Build the mapped tree
    map_tree = get_geo_map_tree(
        geo_ccd_map,
        geo_ann_str="type",
        taxon_map=taxon_map,
        branch_length_map=branch_length_map
    )

    if output_file:
        click.echo("We will write the output file to a tree, WIP, currently nothing happens.")
    else:
        # Print in format 5 with features
        print(map_tree.write(format=5, features=["type"], format_root_node=True))


if __name__ == "__main__":
    tree_file = (f"{Path(__file__).parent.absolute().parent.parent.parent.parent}/examples/"
                 f"data/h3n2-bdmm.h3n2_2deme.trees")

    # todo missing output file...
    main(
        [
            "--trees-file", tree_file,
            "--ccd-type", "ccd0",
            "--burn-in", "0.1"
        ],
    )
