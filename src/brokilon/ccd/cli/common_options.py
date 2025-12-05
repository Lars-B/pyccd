import click


def common_options(func):
    """
    Reusable Click options for ccd CLIs that are shared among all CLIs
    """
    func = click.option(
        "--trees-file",
        required=True,
        type=click.Path(exists=True, dir_okay=False),
        help="Path to the tree file",
    )(func)
    func = click.option(
        "--ccd-type",
        type=click.Choice(("ccd0", "ccd1", "0", "1"), case_sensitive=False),
        default="ccd1",
        callback=validate_ccd_type,
        help="Type of CCD to use for probability calculation",
    )(func)
    func = click.option(
        "--burn-in",
        type=float,
        default=0.1,
        show_default=True,
        required=False,
        help="Fraction of trees to skip from the start."
    )(func)
    func = click.option(
        "--output-file",
        type=click.Path(writable=True, dir_okay=False),
        default=None,
        help="Path to save the CCD-MAP tree"
    )(func)
    func = click.option(
        "--verbose",
        is_flag=True,
        default=False,
    )(func)
    return func


def validate_ccd_type(ctx, param, value):
    match value.lower():
        case "ccd0" | "0":
            return 0
        case "ccd1" | "1":
            return 1
        case _:
            raise click.UsageError(f"Invalid CCD type: {value}")
