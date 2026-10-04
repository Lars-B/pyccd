import networkx as nx


def sranges_map_to_networkx(sranges_map):
    """
    Convert sranges_map: Dict[SRangesClade, Dict[AncestralSplit, support]] to NetworkX.

    Structure: Each clade has MULTIPLE possible splits (with support scores).
    This captures phylogenetic uncertainty at each node.
    """
    G = nx.DiGraph()
    visited_clades = set()

    def clade_to_node_id(clade):
        """Create unique, hashable node ID from SRangesClade."""
        taxa = tuple(sorted(clade.clade)) if clade.clade else ('SPECIAL',)
        range_id = clade.ancestral_range or 'NO_RANGE'
        return f"clade_{hash(taxa) % 10000}_{range_id}"

    def split_to_edge_label(split, support):
        """Create informative label for split edges."""
        anc_taxa = len(split.ancestor.clade) if split.ancestor.clade else 0
        desc_taxa = len(split.descendant.clade) if split.descendant.clade else 0
        return f"sup={support}\n({anc_taxa}/{desc_taxa})"

    def add_clade_recursively(clade):
        """Recursively add a clade and all its possible splits."""
        if clade in visited_clades:
            return

        clade_id = clade_to_node_id(clade)

        # Add clade node
        is_terminal = len(clade.clade) <= 1 if clade.clade else True
        G.add_node(
            clade_id,
            clade=clade,
            is_terminal=is_terminal,
            ancestral_range=clade.ancestral_range,
            num_solutions=len(sranges_map.get(clade, {}))
        )

        # If this clade has resolutions (splits), process them
        if clade in sranges_map:
            splits_dict = sranges_map[clade]

            for split, support in splits_dict.items():
                # Get child clades from this split
                anc_clade = split.ancestor
                desc_clade = split.descendant

                anc_id = clade_to_node_id(anc_clade)
                desc_id = clade_to_node_id(desc_clade)

                # Create intermediate node for the split (optional, for clarity)
                split_id = f"split_{clade_id}_to_{hash(id(split)) % 10000}"
                G.add_node(
                    split_id,
                    split=split,
                    support=support,
                    is_split_node=True
                )

                # Connect clade → split → children
                G.add_edge(clade_id, split_id, weight=support)
                G.add_edge(split_id, anc_id, role='ancestor',
                           taxa=len(anc_clade.clade) if anc_clade.clade else 0)
                G.add_edge(split_id, desc_id, role='descendant',
                           taxa=len(desc_clade.clade) if desc_clade.clade else 0)

                # Recurse on child clades (they may have their own splits)
                add_clade_recursively(anc_clade)
                add_clade_recursively(desc_clade)
        else:
            # Terminal clade - no further resolution
            pass

        visited_clades.add(clade)

    # Process all clades - handles forests automatically
    for clade in sranges_map.keys():
        if clade not in visited_clades:
            add_clade_recursively(clade)

    return G
