import click

from brokilon.core import read_nexus_trees
from brokilon.ccd.domain.sranges import sranges, sranges_ccd_graph_visualizer
from brokilon.ccd.cli.common_options import common_options


@click.command()
@common_options
@click.option(
    "--output-file",
    "html_output",
    type=click.Path(dir_okay=False, writable=True),
    default="sranges-ccd-vis.html",
    show_default=True,
    help="Path for the output HTML visualisation.",
)
@click.option(
    "--reverse-taxon-map",
    is_flag=True,
    help="Use the reverse taxon map for display labels.",
)
def main(trees_file, ccd_type, burn_in, output_file, verbose,
         html_output, reverse_taxon_map):
    """Generate an interactive sRanges Cytoscape HTML visualisation."""

    if ccd_type != 1:
        raise click.ClickException(
            "There is currently only a CCD1 implemented."
        )

    if not 0.0 <= burn_in < 1.0:
        raise click.BadParameter(
            "Must be between 0.0 (inclusive) and 1.0 (exclusive).",
            param_hint="--burn-in",
        )

    if verbose:
        click.echo("Parsing input trees...", err=True)

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
        click.echo(f"After burn-in: {len(trees)} trees.", err=True)

    clade_counts, clade_split_counts, sranges_set, sampled_ancestors = (
        sranges.get_sranges_map(trees, taxon_map)
    )

    reverse_taxon_map_dict = {
        value: key for key, value in taxon_map.items()
    }

    if verbose:
        click.echo("Generating Cytoscape HTML visualisation...", err=True)

    sranges_ccd_graph_visualizer.sranges_map_to_cytoscape_html(
        clade_split_counts,
        reverse_taxon_map_dict,
        clade_counts,
        output_file=html_output,
        use_reverse_taxon_map=reverse_taxon_map,
    )

    if verbose:
        click.echo(f"Visualisation written to {html_output}", err=True)


if __name__ == "__main__":
    main()
