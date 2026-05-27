import csv
import datetime
import datetime as dt
from pathlib import Path

import click
import pandas as pd
from brokilon.ccd.domain.transmission import read_breath_nexus
from brokilon.ccd.domain.transmission.find_infectors import (find_infector_unknown,
                                                             find_infector_with_data,
                                                             find_infector)


class DateExtractionError(Exception):
    """Raised if taxon names do not contain a date"""
    pass


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
            except Exception:
                continue
    raise DateExtractionError(
        f"Could not extract date from label {taxon_label}.\n"
        f"Input for separation was '{sep}' and the date format was '{fmt}'."
    )


def get_root_age_from_leafs(tree, taxon_map, sep, fmt, scale):
    root_node = tree.get_tree_root()
    # scaling floats to dates
    scaling_list = []
    for l in tree:
        cur_root_dist = l.get_distance(root_node)
        cur_date = extract_date_from_label(taxon_map[int(l.name)], sep, fmt)
        scaling_list.append((cur_root_dist, cur_date))

    root_dates = []
    for root_dist, date in scaling_list:
        root_dates.append(date - datetime.timedelta(days=root_dist * scale))

    unique_dates = sorted(set(root_dates))
    if len(unique_dates) > 1:
        diff_days = (unique_dates[-1] - unique_dates[0]).days
        if diff_days > 10:
            print(f"Root date range: {unique_dates[0]} to {unique_dates[-1]}. ({diff_days} days)")
            raise ValueError("This is too much!?")

    # print("Root date that will be used:", unique_dates[0])
    return unique_dates[0]


def get_root_age_with_date(tree, start_date, scale, taxon_map):
    """
    Using the start date on the oldest taxon we can extract the root age using the scale and float
    conversion

    :param tree: A tree to get a root date for
    :param start_date: The date asssumed for the most recent leaf (furthest from root)
    :param scale: Scale for 1.0 float to days/years
    :param taxon_map: The corresponding taxon map
    :return:
    """
    root_node = tree.get_tree_root()
    furthest_root_distance = max(l.get_distance(root_node) for l in tree)
    root_date = start_date - datetime.timedelta(days=furthest_root_distance * scale)

    # Writing a dataframe of dates of all leafs...
    date_list = [
        (
            taxon_map[int(l.name)],
            float_to_date(root_date, l.get_distance(root_node), scale)
        )
        for l in tree
    ]

    return root_date, date_list


def float_to_date(root_date, float_val, scale):
    days_offset = float_val * scale
    return root_date + datetime.timedelta(days=days_offset)


def translate(value, taxon_map):
    if str(value).startswith("Unknown"):
        return str(value)
    if str(value).startswith("block"):
        return str(value)
    try:
        return taxon_map.get(int(value), "Unknown_?")
    except (ValueError, TypeError):
        return "Unknown_?"


def extracting_data(tree, taxon_map, sep, fmt, scale):
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
    # unique_data = [x for x in {tuple(sublist) for sublist in data_frame}]
    unique_data = list({tuple(sublist) for sublist in data_frame})

    scaled_data_frame = []

    # This is for leaf to date map data if leaf labels don't have dates in them
    leaf_dates = None

    try:
        root_date = get_root_age_from_leafs(tree, taxon_map, sep, fmt, scale)
    except DateExtractionError:
        start_leaf_date = datetime.date.today()
        root_date, leaf_dates = get_root_age_with_date(tree, start_leaf_date, scale, taxon_map)

    for infector, infectee, start, blockcount in unique_data:
        scaled_data_frame.append([
            translate(infector, taxon_map),
            translate(infectee, taxon_map),
            float_to_date(root_date, start, scale),
            # todo we can translate blockcount to three types of infection events...
            blockcount if blockcount else "NaN",
        ])

    # print(data_frame)
    # print(scaled_data_frame)
    df = pd.DataFrame(scaled_data_frame, columns=["Infector", "Infectee", "Infection Start",
                                                  "Blockcount"])
    return df, leaf_dates


def process_trees_file(
        trees_file,
        output=None,
        burn_in=0.1,
        date_sep="+",
        date_format="%Y-%m-%d",
        scale=365.24219,
        progress_callback=None,
        log_callback=None,
):
    trees_file = Path(trees_file).absolute()

    if not 0.0 <= burn_in < 1.0:
        raise ValueError("Burn-in must be between 0.0 (inclusive) and 1.0 (exclusive).")

    from brokilon.core.read_nexus import count_trees_in_nexus

    total_trees = count_trees_in_nexus(trees_file)
    burn_in_end = int(burn_in * total_trees)

    if total_trees - burn_in_end == 0:
        raise ValueError("No trees left after burn-in, reduce value of burn_in.")

    if log_callback:
        log_callback(
            f"File contains {total_trees} trees. "
            f"Removing up to tree {burn_in_end} as burnin."
        )

    trees, taxon_map = read_breath_nexus(
        trees_file,
        parse_taxon_map=True,
        burn_in=burn_in
    )

    if log_callback:
        log_callback(f"Parsed {len(trees)} trees.")

    all_results = []
    all_leaf_dates = []

    for i, tree in enumerate(trees):
        if progress_callback:
            progress_callback(i, len(trees))

        cur_df, leaf_dates = extracting_data(
            tree, taxon_map, date_sep, date_format, scale
        )

        cur_df["tree_index"] = i + burn_in_end + 1
        all_results.append(cur_df)
        all_leaf_dates.append(leaf_dates)

    final_df = pd.concat(all_results, ignore_index=True)

    first_leaf_dates = all_leaf_dates[0] if all_leaf_dates else None

    if first_leaf_dates:
        inconsistent = any(ld != first_leaf_dates for ld in all_leaf_dates[1:])
        if inconsistent and log_callback:
            log_callback("Warning: leaf dates differ across trees!")

    # OUTPUT HANDLING (kept simple for GUI reuse)
    if output:
        output_path = Path(output).absolute()
        final_df.to_csv(output_path, index=False)

        leaf_dates_path = output_path.with_name(output_path.stem + "_leaf_dates.csv")

        if first_leaf_dates:
            with open(leaf_dates_path, "w", newline="\n") as f:
                writer = csv.writer(f)
                for i, leafs in enumerate(all_leaf_dates):
                    for taxon, d in leafs:
                        writer.writerow((taxon, d.isoformat(), i))

        return final_df, str(output_path), str(leaf_dates_path)

    return final_df, None, None


@click.command()
@click.option("--trees-file", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--output", type=click.Path(writable=True, dir_okay=False), default=None)
@click.option("--burn-in", type=float, default=0.1, show_default=True)
@click.option("--date-sep", type=str, default="+", show_default=True)
@click.option("--date-format", type=str, default="%Y-%m-%d", show_default=True)
@click.option("--scale", type=float, default=365.24219, show_default=True)
def main(trees_file, output, burn_in, date_sep, date_format, scale):
    def log(msg):
        click.echo(msg, err=True)

    def progress(i, total):
        # CLI-friendly minimal progress
        if i % 10 == 0:
            click.echo(f"Processing {i}/{total}", err=True)

    process_trees_file(
        trees_file=trees_file,
        output=output,
        burn_in=burn_in,
        date_sep=date_sep,
        date_format=date_format,
        scale=scale,
        log_callback=log,
        progress_callback=progress,
    )


if __name__ == '__main__':
    main(
        args=["--trees-file", "../../../../../../testing/truth.trees",
              "--burn-in", "0",
              "--output", "../../../../../../testing/refactor_out.log"],
    )
