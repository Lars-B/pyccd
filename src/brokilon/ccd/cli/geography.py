import click
from pathlib import Path
from brokilon.core.read_nexus import read_nexus_trees
from brokilon.ccd.domain.phylogeography import get_geo_map, get_geo_map_tree
from brokilon.ccd.cli.common_options import common_options


@click.command()
@common_options
def main(trees_file, ccd_type, burn_in, output_file, verbose):
    """
    Process a Nexus tree file and generate a geo-mapped tree.
    """
    trees_file = Path(trees_file).absolute()

    if verbose:
        click.echo("Reading Nexus trees...")

    # todo optimize the reading of the file, takes a looooot of time atm...

    # Read trees and taxon map
    trees, taxon_map = read_nexus_trees(trees_file, parse_taxon_map=True)

    # Skip a fraction of trees if needed
    start_idx = int(len(trees) * burn_in)
    trees = trees[start_idx:]

    if verbose:
        click.echo("Computing CCD graph and probabilities...")

    # todo this needs to be optimized or somehow made to work, shouldn't take thiiiis long?
    # Compute geo map and branch lengths
    geo_ccd_map, branch_length_map, clade_count_map = get_geo_map(
        trees, geo_ann_str="type", ccd_type=ccd_type
    )

    # todo look into ccd0 for this and see if I can optimize it a bit more...

    if verbose:
        click.echo("Computing CCD-MAP tree...")

    # Build the mapped tree
    map_tree = get_geo_map_tree(
        geo_ccd_map,
        geo_ann_str="type",
        taxon_map=None if output_file else taxon_map,
        branch_length_map=branch_length_map,
        clade_count_map=clade_count_map
    )

    # todo possibly rename support to ccd-prob or something...
    newick_map = map_tree.write(format=5, features=["type", "support"], format_root_node=True)

    # todo fix naming here, output file and output file stream, check for if exists stuff etc...
    if output_file:
        click.echo("Writing MAP tree to file...")

        with open(output_file, "x") as output_file_stream:
            with open(trees_file, "r") as input_trees_file:
                for line in input_trees_file:
                    if line.strip().startswith("tree "):
                        break
                    output_file_stream.write(line)

                output_file_stream.write(f"tree ext_CCD_MAP = {newick_map}\nEnd;\n")
    else:
        # Print in format 5 with features
        # print(map_tree.write(format=5, features=["type"], format_root_node=True))
        click.echo(newick_map)


if __name__ == "__main__":
    tree_file = (f"{Path(__file__).parent.absolute().parent.parent.parent.parent}/examples/"
                 f"data/h3n2-bdmm.h3n2_2deme.trees")
    out_file = (f"{Path(__file__).parent.absolute().parent.parent.parent.parent}/examples/"
                 f"data/ext_ccd_summary.tree")

    main(
        [
            "--trees-file", tree_file,
            "--output-file", out_file,
            "--ccd-type", "ccd0",
            "--burn-in", "0.1"
        ],
    )
