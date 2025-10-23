import re


def change_taxon_map(input_nexus, output_nexus, new_map):
    # Function that will apply the new map to the input nexus and save it as output_nexus

    # todo check file exists
    # todo check that the map is compatible etc...

    begin_map = re.compile('\t?translate\n', re.I)
    end = re.compile('\t*?;\n?')
    re_tree = re.compile("\t?tree .*=? (.*$)", flags=re.I | re.MULTILINE)
    re_taxa = re.compile('([0-9]+)([\\[:])')

    new_map_reversed = {v: k for (k, v) in new_map.items()}
    old_map = {}
    within_map = False
    finished_map = False

    with open(input_nexus, "r") as in_file, open(output_nexus, "w+") as out_file:
        for line in in_file:
            if begin_map.match(line):
                # we enter the taxon map lines
                out_file.write(line)
                within_map = True
            if within_map and end.match(line):
                # The taxon map is finished
                out_file.write(";\n")
                within_map = False
                finished_map = True
            if within_map and not begin_map.match(line):
                # We are within the taxon map, all the line need to be extracted
                split = line.split()
                cur_key = int(split[0])
                old_map[cur_key] = split[1][:-1] if split[1][-1] == "," else split[1]
                out_file.write(f"\t\t{cur_key} {new_map[cur_key]},\n")
            if not finished_map and not within_map:
                # Write everything that comes before the taxon map to the new file
                out_file.write(line)
            if re_tree.match(line):
                # matching a tree, need to change the taxon integer matches accordingly
                # extracting the tree newick string
                tree_string = re.split(re_tree, line)[1]
                # apply new taxon map ...
                new_newick = re_taxa.sub(lambda m: m.group().replace(m.group(1),
                                                                     str(new_map_reversed[
                                                                             old_map[int(m.group(1))]])),
                                         tree_string)
                out_file.write(f"{line.split('=')[0]}= {new_newick}\n")
        out_file.write("End;")


def get_mapping_dict(file: str) -> dict:
    """
    Returns the taxon mapping of the nexus file as a dictionary

    :param file: A nexus file path
    :type file: str
    :return: Dictionary containing the mapping of taxa(values) to int(keys)
    :rtype: dict {int --> str}
    """

    begin_map = re.compile('\t?translate\n', re.I)
    end = re.compile('\t*?;\n?')

    mapping = {}

    begin = False
    with open(file) as f:
        for line in f:
            if begin:
                if end.match(line):
                    break
                split = line.split()
                mapping[int(split[0])] = split[1][:-1] if split[1][-1] == "," else split[1]

            if begin_map.match(line):
                begin = True
    return mapping
