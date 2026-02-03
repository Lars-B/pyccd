import click

from brokilon.ccd.domain.transmission import read_breath_nexus
from brokilon.core import read_nexus_trees
from brokilon.majority_rule_consensus.majority_consensus import (regular_mrc,
                                                                 phygeo_mrc, transmission_mrc)


@click.group()
@click.option(
    "--trees-file",
    "trees_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the tree file",
)
@click.option(
    "--burn-in",
    "burn_in",
    type=float,
    default=0.1,
    show_default=True,
    required=False,
    help="Fraction of trees to skip from the start."
)
@click.option(
    "--major-thr",
    "major_thr",
    type=float,
    default=0.5,
    show_default=True,
    required=False,
    help="Fraction of trees that should contain the clades."
)
@click.option(
    "--output-file",
    "output_file",
    type=click.Path(writable=True, dir_okay=False),
    default=None,
    help="Path to save the MRC tree."
)
@click.pass_context
def mrc(ctx, trees_file, burn_in, major_thr, output_file):
    """Compute majority rule consensus tree."""
    ctx.ensure_object(dict)
    ctx.obj["trees_file"] = trees_file
    ctx.obj["burn_in"] = burn_in
    ctx.obj["major_thr"] = major_thr
    ctx.obj["output_file"] = output_file


@mrc.command()
@click.pass_context
def regular(ctx):
    """Regular majority rule consensus (default)."""
    trees_file = ctx.obj["trees_file"]
    burnin = ctx.obj["burn_in"]
    major_thr = ctx.obj["major_thr"]
    click.echo(f"Regular MRC on {trees_file} with burnin {burnin}")
    trees, map = read_nexus_trees(trees_file, burn_in=burnin, parse_taxon_map=True)
    mrc_tree = regular_mrc(trees)
    # todo finish writing to console, will need file output argumetn too


@mrc.command()
@click.option(
    "--annotation",
    default="type",
    show_default=True,
    help="Node annotation to use for geographic MRC.",
)
@click.pass_context
def geo(ctx, annotation):
    """Geographic MRC (supports --annotation option)."""
    trees_file = ctx.obj["trees_file"]
    burnin = ctx.obj["burn_in"]
    major_thr = ctx.obj["major_thr"]
    output_file = ctx.obj["output_file"]
    click.echo(
        f"Geo MRC on {trees_file} with burnin {burnin}, annotation={annotation}"
    )

    # todo geo
    trees, map = read_nexus_trees(trees_file, burn_in=burnin, parse_taxon_map=True)
    mrc_tree = phygeo_mrc(trees, annotation, major_thr=major_thr)

    newick_mrc = mrc_tree.write(format=5, features=["type", "support"], format_root_node=True)

    if output_file:
        click.echo("Writing MRC tree to file...")

        with open(output_file, "w") as out_tfile:
            with open(trees_file, "r") as in_tfile:
                for line in in_tfile:
                    if line.strip().startswith("tree "):
                        break
                    out_tfile.write(line)

            out_tfile.write(f"tree MRC = {newick_mrc}\nEnd;\n")
    else:
        click.echo(newick_mrc)



@mrc.command(name="transm")
@click.pass_context
def transmission(ctx):
    """Transmission MRC."""
    trees_file = ctx.obj["trees_file"]
    burnin = ctx.obj["burn_in"]
    major_thr = ctx.obj["major_thr"]
    click.echo(f"Transmission MRC on {trees_file} with burnin {burnin}")

    # todo tranmission
    trees, map = read_breath_nexus(trees_file, burn_in=burnin, parse_taxon_map=True)
    mrc_tree = transmission_mrc(trees)

    # todo output


if __name__ == '__main__':
    from pathlib import Path

    tree_file = (f"{Path(__file__).parent.absolute().parent.parent.parent}/examples/"
                 f"data/h3n2-bdmm.h3n2_2deme.trees")
    # out_file = (f"{Path(__file__).parent.absolute().parent.parent.parent.parent}/examples/"
    #             f"data/ext_ccd_summary.tree")

    mrc(
        [
            "--trees-file", tree_file,
            "geo",
            "--annotation", "type",
        ],
    )
