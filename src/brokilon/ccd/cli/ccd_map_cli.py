import click

from brokilon.core import read_nexus_trees


@click.command()
@click.option(
    "--trees-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the tree file",
)
@click.option(
    "--output-file",
    type=click.Path(writable=True, dir_okay=False),
    default=None,
    help="Path to save the CCD-MAP tree"
)
@click.option(
    "--burn-in",
    type=float,
    default=0.1,
    required=False,
    is_eager=True,
    show_default=True,
    help="Burn-in proportion between 0.0 and 1.0."
)
@click.option(
    "--ccd-type",
    type=click.Choice(("ccd1", "ccd0", "0", "1"), case_sensitive=False),
    required=False,
    default="ccd1",
    help="Something...."
)
def main(trees_file, output_file, burn_in, ccd_type):

    trees = read_nexus_trees(trees_file, breath_trees=False, label_transm_history=False)
    trees = trees[int(len(trees)*burn_in):]

    if len(trees) == 0:
        click.UsageError(f"There are no trees left after burn-in!"
                         f" Try reducing it... (currently = {burn_in})")

    click.echo(f"Parsed {len(trees)} trees.", err=True)

    # todo check which type of ccd we have...
    match ccd_type:
        case "ccd0" | "0":
            click.echo("Type is ccd0")
        case "ccd1" | "1":
            click.echo("Default option is ccd1")
        case _:
            raise click.UsageError(f"Invalid CCD type: {ccd_type}")

    if output_file:
        # write to file
        click.echo("Writing file")
    else:
        # print to sys.stdout
        click.echo("Writing to std:out")

    # todo figure out annotation process of the tree
    return 0


if __name__ == '__main__':
    main()
