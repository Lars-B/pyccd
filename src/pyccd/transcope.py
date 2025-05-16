import os
import argparse
import sys

from src.pyccd.read_nexus import read_nexus_trees
from src.pyccd.transmission_ccd import get_transmission_maps, get_transmission_ccd_tree_bottom_up


def main():
    parser = argparse.ArgumentParser(description='Transmission CCD MAP tree computation',
                                     prog='transcope')
    parser.add_argument(
        '-i', '--input-trees',
        help='Input tree file path',
        required=True,
    )
    parser.add_argument('-o', '--output-tree', nargs='?',
                        help='Output Tree file (default: %(default)s)',
                        type=argparse.FileType('w'), default=sys.stdout,
                        )
    parser.add_argument('-b', '--burn-in', nargs='?',
                        default=0.0,  # default value
                        const=0.1,  # default value when just --burn-in is passed
                        type=float,
                        help='Burn-in proportion between 0.0 and 1.0 (default: %(default)s))'
                        )
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    input_trees = args.input_trees
    output_tree = args.output_tree
    burn_in = args.burn_in
    overwrite = args.overwrite
    verbose = args.verbose

    if not os.path.isfile(input_trees):
        print(f"Error: input tree file {input_trees} does not exist!", file=sys.stderr)
        sys.exit(1)
    if os.path.exists(output_tree.name) and not overwrite:
        print("Output tree already exists. Use '--overwrite' to overwrite the file.",
              file=sys.stderr)
        sys.exit(1)
    if not 0.0 <= burn_in < 1.0:
        print("Burn-in must be between 0.0 (inclusive) and 1.0 (exclusive).", file=sys.stderr)
        sys.exit(1)

    trees = read_nexus_trees(input_trees, breath_trees=True)
    trees = trees[int(burn_in * len(trees)):]
    if len(trees) < 1:
        print("Input trees empty after burn-in removal, maybe burn-in too high?", file=sys.stderr)
        sys.exit(1)

    m1, m2, blockcount_map, branch_lengths_map = get_transmission_maps(trees)
    newick_map = get_transmission_ccd_tree_bottom_up(m1, m2, blockcount_map, branch_lengths_map)

    with open(input_trees, 'r', encoding="UTF-8") as infile:
        for line in infile:
            if line.strip().startswith("tree "):
                break  # Stop reading when the tree section starts
            output_tree.write(line)

        output_tree.write(f"tree tCCD_MAP = {newick_map};\nEnd;\n")

    # todo verbose needs to be used with progress etc...
    # todo empty run of the program should either default ot print help or print a meaningful error?
    print("Done invoking transcope.")
    return 0


if __name__ == '__main__':
    main()
