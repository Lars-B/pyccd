from brokilon.core import read_nexus_trees
from brokilon.ccd.domain.transmission.label_transmission_history import label_transmission_tree


def read_breath_nexus(file, parse_taxon_map: bool = False, burn_in: float = 0) -> tuple:
    results = read_nexus_trees(file, parse_taxon_map=parse_taxon_map, burn_in=burn_in)

    if parse_taxon_map:
        trees, taxon_map = results
    else:
        trees = results
        taxon_map = None

    for tree in trees:
        label_transmission_tree(tree)

    if parse_taxon_map:
        return trees, taxon_map
    return trees
