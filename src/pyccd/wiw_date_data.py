import sys
from datetime import timedelta

import click
import pandas as pd
import datetime as dt

from read_nexus import read_nexus_trees
from pyccd.find_infectors import find_infector_unknown, find_infector_with_data, find_infector

global SCALE
# Days per year, i.e. 1.0 float of branch length equals this value
SCALE = 365.24219


def extract_date_from_label(taxon_label: str,
                            sep: str = "+",
                            fmt: str = "%Y-%m-%d"):
    try:
        date_str = taxon_label.split(sep)[1]
        return dt.datetime.strptime(date_str, fmt).date()
    except Exception:
        pass

    fallback_separators = ("+", ":", "-", "_",)
    fallback_formats = ("%Y-%m-%d", "%d-%m-%Y", "%Y/%m/%d", "%m/%d/%Y",)

    for s in fallback_separators:
        if s not in taxon_label:
            continue
        parts = taxon_label.split(s)
        if len(parts) < 2:
            continue
        for f in fallback_formats:
            try:
                return dt.datetime.strptime(parts[1], f).date()
            except Exception as e:
                continue
    raise ValueError(f"Could not extract date from label {taxon_label}.\n"
                     f"Input for separation was '{sep}' and the date format was '{fmt}'.")


def get_root_age_from_leafs(tree, taxon_map, sep, fmt):
    root_node = tree.get_tree_root()
    # scaling floats to dates
    scaling_list = []
    for l in tree:
        cur_root_dist = l.get_distance(root_node)
        cur_date = extract_date_from_label(taxon_map[int(l.name)], sep, fmt)
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
        return str(value)
    if str(value).startswith("block"):
        return str(value)
    try:
        return taxon_map.get(int(value), "Unknown_?")
    except (ValueError, TypeError):
        return "Unknown_?"


def extracting_data(tree, taxon_map, sep, fmt):
    root_node = tree.get_tree_root()

    data_frame = []
    unknown_nodes_todo = []
    seen_unknown_labels = set()

    for leaf in tree:
        og_infector = find_infector(leaf)
        infector, cur_data, unknown_node = find_infector_with_data(leaf, root_node)

        if unknown_node is not None and unknown_node.transm_ancest not in seen_unknown_labels:
            unknown_nodes_todo.append(unknown_node)
            seen_unknown_labels.add(unknown_node.transm_ancest)

        data_frame.extend(cur_data)
        if og_infector != infector:
            raise ValueError("This should not happen?!")

    for todo_node in unknown_nodes_todo:
        assert todo_node.transm_ancest.startswith("Unknown"), "Something went wrong!"
        cur_data = find_infector_unknown(todo_node, root_node)
        data_frame.extend(cur_data)

    # print("-----")
    # for i in data_frame:
    #     print(i)
    # print("-----")
    # todo not sure why some nodes are not unique but we can just remove the duplicate events....
    unique_data = [x for x in {tuple(sublist) for sublist in data_frame}]


    scaled_data_frame = []
    root_date = get_root_age_from_leafs(tree, taxon_map, sep, fmt)
    for infector, infectee, start, blockcount in unique_data:
        scaled_data_frame.append([
            translate(infector, taxon_map),
            translate(infectee, taxon_map),
            float_to_date(root_date, start),
            # todo we can translate blockcount to three types of infection events...
            blockcount if blockcount else "NaN",
        ])

    # print(data_frame)
    # print(scaled_data_frame)
    df = pd.DataFrame(scaled_data_frame, columns=["Infector", "Infectee", "Infection Start",
                                                  "Blockcount"])
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
@click.option(
    "--burn-in",
    type=float,
    default=0.1,       # default if option is NOT passed
    required=False,
    is_eager=True,
    show_default=True,
    help="Burn-in proportion between 0.0 and 1.0."
)
@click.option(
    "--date-sep",
    type=str,
    default="+",
    show_default=True,
    help="Separator used in taxon labels to split ID and date."
)
@click.option(
    "--date-format",
    type=str,
    default="%Y-%m-%d",
    show_default=True,
    help="Date format string to parse dates."
)
def main(trees_file, output, burn_in, date_sep, date_format):
    trees, taxon_map = read_nexus_trees(
        trees_file,
        breath_trees=True,
        parse_taxon_map=True
    )

    burnin_index = int(burn_in * len(trees))
    trees = trees[burnin_index:]

    click.echo(f"Parsed {len(trees)} trees.", err=True)

    all_results = []

    with (click.progressbar(enumerate(trees, start=burnin_index),
                            length=len(trees),
                            label="Processing trees") as bar):
        for i, tree in bar:
            cur_df = extracting_data(tree, taxon_map, date_sep, date_format)
            cur_df["tree_index"] = i  # add the tree index as a new column
            all_results.append(cur_df)

    final_df = pd.concat(all_results, ignore_index=True)

    if output:
        final_df.to_csv(output, index=False)
        click.echo(f"Saved results to {output}", err=True)
    else:
        final_df.to_csv(sys.stdout, index=False)


if __name__ == '__main__':
    main(
        # args=["--trees-file", "../../tests/data/Testing_date_extraction.trees",
        args=["--trees-file", "../../tests/data/small_example.trees",
              "--burn-in", "0.0"],
    )
