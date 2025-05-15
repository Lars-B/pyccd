import os
import argparse


def main():
    parser = argparse.ArgumentParser(description='Transmission CCD MAP CLI')
    parser.add_argument('-i', '--input-trees', nargs='+',)
    parser.add_argument('-o', '--output-tree', nargs='+',)

    args = parser.parse_args()
    input_trees = args.input_trees
    output_tree = args.output_tree
    print(input_trees)
    print(output_tree)
    # todo parsing all the arguments
    print("---")
    return 0


if __name__ == '__main__':
    main()
