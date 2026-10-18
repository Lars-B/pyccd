"""Generate an interactive Cytoscape HTML visualization of the stratigraphic-range CCD graph."""

from pathlib import Path

import click

from brokilon.ccd.domain.sranges import sranges, sranges_ccd_graph_visualizer
from brokilon.core import read_nexus_trees


@click.command()
@click.option(
    "--trees-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the input Nexus tree file.",
)
@click.option(
    "--burn-in",
    type=click.FloatRange(min=0.0, max=1.0, max_open=True),
    default=0.1,
    show_default=True,
    help="Fraction of initial trees to discard.",
)
@click.option(
    "--output-file",
    type=click.Path(dir_okay=False, writable=True),
    default="sranges-ccd-vis.html",
    show_default=True,
    help="Path for the output HTML visualization.",
)
@click.option(
    "--reverse-taxon-map/--no-reverse-taxon-map",
    default=False,
    show_default=True,
    help="Use the reverse taxon map to display only integers,"
         " without you will get full taxon labels.",
)
@click.option(
    "--verbose",
    is_flag=True,
    help="Print progress messages.",
)
def main(trees_file, burn_in, output_file, reverse_taxon_map, verbose):
    """Generate an interactive fossil stratigraphic-range visualization."""

    if output_file and Path(output_file).exists():
        raise click.ClickException(
            f"Output file already exists: {output_file}"
        )

    if verbose:
        click.echo("Parsing input trees...", err=True)

    trees, taxon_map = read_nexus_trees(
        trees_file,
        parse_taxon_map=True,
    )

    # Remove burn-in.
    trees = trees[int(burn_in * len(trees)):]

    if not trees:
        raise click.ClickException(
            "Input trees empty after burn-in removal. "
            "Consider lowering --burn-in."
        )

    if verbose:
        click.echo(f"After burn-in: {len(trees)} trees.", err=True)

    clade_counts, clade_split_counts, sranges_set, sampled_ancestors = (
        sranges.get_sranges_map(trees, taxon_map)
    )

    reverse_taxon_map_dict = {
        value: key for key, value in taxon_map.items()
    }

    if verbose:
        click.echo("Generating HTML visualization...", err=True)

    sranges_ccd_graph_visualizer.sranges_map_to_cytoscape_html(
        clade_split_counts,
        reverse_taxon_map_dict,
        clade_counts,
        output_file=output_file,
        use_reverse_taxon_map=reverse_taxon_map,
    )

    if verbose:
        click.echo(f"Visualization written to {output_file}", err=True)


if __name__ == "__main__":
    main()
