"""
Module contains the transcope Command line interface
"""
import sys
from pathlib import Path

import click

from brokilon.ccd.cli.common_options import common_options
from brokilon.ccd.domain.transmission import (get_transmission_maps,
                                              get_transmission_ccd_tree_bottom_up)
from brokilon.ccd.domain.transmission import read_breath_nexus


@click.command()
@common_options
@click.option(
    "--seed",
    type=int,
    default=1337,
    help="Random number seed for CCD-MAP tree tiebreaking"
)
def main(trees_file, ccd_type, burn_in, output_file, verbose, seed):
    """
    Command line interface to calculate a transmission CCD-MAP tree with input options.
    """

    trees_file = Path(trees_file).absolute()

    if not 0.0 <= burn_in < 1.0:
        print("Burn-in must be between 0.0 (inclusive) and 1.0 (exclusive).", file=sys.stderr)
        sys.exit(1)

    if verbose:
        print("Parsing input trees...", file=sys.stderr)

    trees, taxon_map = read_breath_nexus(trees_file, True)
    trees = trees[int(burn_in * len(trees)):]

    if len(trees) < 1:
        print("Input trees empty after burn-in removal, maybe burn-in too high?", file=sys.stderr)
        sys.exit(1)
    if verbose:
        print(f"After burn-in there are {len(trees)} trees left...", file=sys.stderr)

    if ccd_type == 1:
        m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(
            trees,
            type_str="Ancestry"
        )
    elif ccd_type == 0:
        raise NotImplementedError("Currently not implemented in this CLI")

    # todo add similar to geography option to make just a newick string that will be printed
    newick_map = get_transmission_ccd_tree_bottom_up(
        m1, m2,
        blockcount_map, branch_lengths_map,
        seed=seed
    )

    if verbose:
        print("Writing transmission CCD-MAP tree to file...", file=sys.stderr)

    if output_file:
        with open(output_file, "x") as output_file_stream:
            with open(trees_file, 'r', encoding="UTF-8") as infile:
                for line in infile:
                    if line.strip().startswith("tree "):
                        break  # Stop reading when the tree section starts
                    output_file_stream.write(line)

                output_file_stream.write(f"tree tCCD_MAP = {newick_map};\nEnd;\n")
    else:
        raise NotImplementedError("No outputfile provided, unsupported atm.")

    if verbose:
        print("Done invoking transcope.", file=sys.stderr)
    return 0


if __name__ == '__main__':
    in_trees = (f"{Path(__file__).parent.absolute().parent.parent.parent.parent}/examples/"
                f"data/breath32sim.trees")

    out_tree = (f"{Path(__file__).parent.absolute().parent.parent.parent.parent}/examples/"
                f"data/tccd.tree")

    main(
        [
            "--trees-file", in_trees,
            "--output-file", out_tree,
            "--verbose",
            "--burn-in", 0.1
        ]
    )
