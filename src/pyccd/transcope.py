"""
Module contains the transcope Command line interface
"""
import os
import argparse
import sys

from .read_nexus import read_nexus_trees
from .transmission_ccd import get_transmission_maps, get_transmission_ccd_tree_bottom_up


def main():
    """
    Command line interface to calculate a transmission CCD-MAP tree with input options.
    """
    prelim_parser = argparse.ArgumentParser(add_help=False)
    prelim_parser.add_argument('--overwrite', action='store_true',)
    prelim_args, _ = prelim_parser.parse_known_args()

    parser = argparse.ArgumentParser(description='Transmission CCD MAP tree computation',
                                     prog='transcope')
    parser.add_argument(
        '-i', '--input-trees',
        help='Input tree file path',
        required=True,
    )
    parser.add_argument('-o', '--output-tree',
                        help='Output Tree file (default: standard output)',
                        type=argparse.FileType(f"{'w' if prelim_args.overwrite else 'x'}"),
                        default=sys.stdout,
                        )
    parser.add_argument('-b', '--burn-in', nargs='?',
                        default=0.0,  # default value
                        const=0.1,  # default value when just --burn-in is passed
                        type=float,
                        help='Burn-in proportion between 0.0 and 1.0 (default: %(default)s))'
                        )
    parser.add_argument('-v', '--verbose',
                        action='store_true', help="Enable verbose status output")
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing output file')
    args = parser.parse_args()

    if not os.path.isfile(args.input_trees):
        print(f"Error: input tree file {args.input_trees} does not exist!", file=sys.stderr)
        sys.exit(1)
    if not 0.0 <= args.burn_in < 1.0:
        print("Burn-in must be between 0.0 (inclusive) and 1.0 (exclusive).", file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print("Parsing input trees...", file=sys.stderr)
    trees = read_nexus_trees(args.input_trees, breath_trees=True)
    trees = trees[int(args.burn_in * len(trees)):]
    if len(trees) < 1:
        print("Input trees empty after burn-in removal, maybe burn-in too high?", file=sys.stderr)
        sys.exit(1)
    if args.verbose:
        print(f"After burn-in there are {len(trees)} trees left...", file=sys.stderr)

    # todo there is no way of controlling which type of tCCD this constructs...
    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees)
    newick_map = get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)

    if args.verbose:
        print("Writing transmission CCD-MAP tree to file...", file=sys.stderr)

    with open(args.input_trees, 'r', encoding="UTF-8") as infile:
        for line in infile:
            if line.strip().startswith("tree "):
                break  # Stop reading when the tree section starts
            args.output_tree.write(line)

        args.output_tree.write(f"tree tCCD_MAP = {newick_map};\nEnd;\n")

    if args.verbose:
        print("Done invoking transcope.", file=sys.stderr)
    return 0


if __name__ == '__main__':
    main()
