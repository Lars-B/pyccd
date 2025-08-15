import sys
from datetime import timedelta

import click
import pandas as pd

from . import read_nexus_trees
from .wiw_network import find_infector_unknown, find_infector_with_data, find_infector

global SCALE
# Days per year, i.e. 1.0 float of branch length equals this value
SCALE = 365.24219


def get_root_age_from_leafs(tree, taxon_map):
    root_node = tree.get_tree_root()
    # scaling floats to dates
    DATE_SEPARATOR = "+"
    DATE_FORMAT = "%Y-%m-%d"
    scaling_list = []
    for l in tree:
        cur_root_dist = l.get_distance(root_node)
        cur_date_str = taxon_map[int(l.name)].split(DATE_SEPARATOR)[1]
        import datetime as dt
        cur_date = dt.datetime.strptime(cur_date_str, DATE_FORMAT).date()
        scaling_list.append((cur_root_dist, cur_date))

    # todo compute scaling from map of root distances to dates
    # scale is in years, days or whatever, can assume 1.0 = 1 year for now
    # from itertools import combinations
    # scales = []
    # for (f1, d1), (f2, d2) in combinations(scaling_list, 2):
    #     float_diff = abs(f2 - f1)
    #     day_diff = abs((d2 - d1).days)
    #
    #     if float_diff != 0:
    #         cur_scale = day_diff / float_diff
    #         scales.append(cur_scale)
    #     else:
    #         print("Problem: Floatdiff is 0.0")
    # avg_scales = sum(scales) / len(scales)
    # variation = max(scales) - min(scales)
    # print(f"The average scale is {avg_scales}")
    # print(f"The variation is {variation}")

    root_dates = []
    for root_dist, date in scaling_list:
        root_dates.append(date - timedelta(days=(root_dist * SCALE)))

    unique_dates = sorted(set(root_dates))
    if len(unique_dates) > 1:
        diff_days = (unique_dates[-1] - unique_dates[0]).days
        if diff_days > 10:
            print(f"Root date range: {unique_dates[0]} to {unique_dates[-1]}. ({diff_days} days)")
            raise ValueError("This is too much!?")

    # print("Root date that will be used:", unique_dates[0])
    return unique_dates[0]


def float_to_date(root_date, float_val):
    days_offset = float_val * SCALE
    return root_date + timedelta(days=days_offset)


def translate(value, taxon_map):
    if str(value).startswith("Unknown"):
        return "Unknown"
    try:
        return taxon_map.get(int(value), "Unknown")
    except (ValueError, TypeError):
        return "Unknown"


def extracting_data(tree, taxon_map):
    root_node = tree.get_tree_root()

    data_frame = []
    unknown_nodes_todo = []

    for leaf in tree:
        og_infector = find_infector(leaf)
        infector, cur_data, unknown_node = find_infector_with_data(leaf, root_node)
        unknown_nodes_todo.append(unknown_node) if unknown_node is not None else None
        data_frame.append(cur_data)
        if og_infector != infector:
            raise ValueError("This should not happen?!")

    for todo_node in unknown_nodes_todo:
        assert todo_node.transm_ancest.startswith("Unknown"), "Something went wrong!"
        cur_data = find_infector_unknown(todo_node, root_node)
        data_frame.extend(cur_data)

    scaled_data_frame = []
    root_date = get_root_age_from_leafs(tree, taxon_map)
    for infector, infectee, start, end, blockcount in data_frame:
        scaled_data_frame.append([
            translate(infector, taxon_map),
            translate(infectee, taxon_map),
            float_to_date(root_date, start),
            float_to_date(root_date, end),
            # todo we can translate blockcount to three types of infection events...
            blockcount
        ])

    # print(data_frame)
    # print(scaled_data_frame)
    df = pd.DataFrame(scaled_data_frame, columns=["Infector", "Infectee", "Infection Start",
                                                  "Infection End", "Infection Type (Blockcount)"])
    return df


@click.command()
@click.option(
    "--trees-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the file containing the trees."
)
@click.option(
    "--output",
    type=click.Path(writable=True, dir_okay=False),
    default=None,
    help="Path to save the CSV. Defaults to stdout."
)
def main(trees_file, output):
    trees, taxon_map = read_nexus_trees(
        trees_file,
        breath_trees=True,
        parse_taxon_map=True
    )
    # removing the first tree
    trees = trees[1:]
    click.echo(f"Parsed {len(trees)} trees.", err=True)

    all_results = []
    # for i, tree in enumerate(trees):
    #     click.echo(f"Processing tree index {i}...", err=True)
    #     cur_df = extracting_data(tree, taxon_map)
    #     all_results.append(cur_df)
    with click.progressbar(trees, label="Processing trees") as bar:
        for tree in bar:
            cur_df = extracting_data(tree, taxon_map)
            all_results.append(cur_df)

    final_df = pd.concat(all_results, ignore_index=True)

    if output:
        final_df.to_csv(output, index=False)
        click.echo(f"Saved results to {output}", err=True)
    else:
        final_df.to_csv(sys.stdout, index=False)


if __name__ == '__main__':
    main()
