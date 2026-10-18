"""Compute the stratigraphic-range CCD-MAP tree from a collection of phylogenetic trees."""

import sys

import click

from brokilon.ccd.cli.common_options import common_options
from brokilon.ccd.domain.sranges import sranges
from brokilon.core import read_nexus_trees


@click.command()
@common_options
@click.option(
    "--debug",
    is_flag=True,
    help="Enable debugging output and tree modifications.",
)
def main(trees_file, ccd_type, burn_in, output_file, verbose, debug):
    """Calculate a stratigraphic-range MAP tree from input Nexus trees."""

    if ccd_type != 1:
        raise NotImplementedError("There is currently only a CCD1 implemented.")

    if not 0.0 <= burn_in < 1.0:
        raise click.BadParameter(
            "Must be between 0.0 (inclusive) and 1.0 (exclusive).",
            param_hint="--burn-in",
        )

    if output_file is None:
        raise click.UsageError("--output-file is required.")

    sranges.DEBUG_TREE = debug

    if verbose:
        print("Parsing input trees...", file=sys.stderr)

    trees, taxon_map = read_nexus_trees(
        trees_file,
        parse_taxon_map=True,
    )

    trees = trees[int(burn_in * len(trees)):]

    if not trees:
        raise click.ClickException(
            "Input trees empty after burn-in removal."
        )

    if verbose:
        print(f"After burn-in: {len(trees)} trees.", file=sys.stderr)

    clade_counts, clade_split_counts, sranges_set, sampled_ancestors = (
        sranges.get_sranges_map(trees, taxon_map)
    )

    seen_resolved, map_tree = sranges.get_sranges_map_tree(
        clade_counts,
        clade_split_counts,
        sranges_set,
        sampled_ancestors,
    )

    nwk_map = map_tree.write(
        format=1,
        format_root_node=True,
        features=["orientation", "ancestral_range"],
        icytree_order=True,
    )

    if verbose:
        print("Writing stratigraphic-range MAP tree...", file=sys.stderr)

    with open(output_file, "x", encoding="UTF-8") as output:
        with open(trees_file, "r", encoding="UTF-8") as infile:
            for line in infile:
                if line.strip().startswith("tree "):
                    break
                output.write(line)

        output.write(f"tree sranges_MAP = {nwk_map};\nEnd;\n")

    if verbose:
        print("Done invoking stratigraphic-range CLI.", file=sys.stderr)


if __name__ == "__main__":
    main()
