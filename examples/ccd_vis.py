import json


def sranges_map_to_cytoscape_html(sranges_map, probability_map=None,
                                  output_file='cytoscape_tree.html'):
    """
    sranges_map structure:
      - Keys: SRangesClade objects
      - Values: Dict[AncestralSplit, support_int]
      - AncestralSplit.ancestor → SRangesClade
      - AncestralSplit.descendant → SRangesClade
    """

    visited_clades = set()
    nodes = []
    edges = []
    seen_edges = set()  # Prevent duplicate edges

    def clade_to_node_id(clade):
        """Create unique, hashable node ID from SRangesClade."""
        if clade.clade:
            taxa = tuple(sorted(clade.clade))
        else:
            taxa = ('SPECIAL_NODE',)
        range_id = clade.ancestral_range or 'NO_RANGE'
        node_id = f"clade_{hash(taxa) % 10000}_{range_id}"
        return node_id

    def add_clade_recursive(clade):
        if clade in visited_clades:
            return

        visited_clades.add(clade)
        clade_id = clade_to_node_id(clade)

        # Determine if this is a terminal clade (leaf)
        is_terminal = len(clade.clade) == 1 if clade.clade else False

        # Add node for this clade
        nodes.append({
            'data': {
                'id': clade_id,
                'label': f"{len(clade.clade) if clade.clade else 0}",
                'terminal': is_terminal,
                'range': clade.ancestral_range or '',
                # 'num_resolutions': len(sranges_map.get(clade, {}))
            }
        })

        # If this clade has resolutions (splits), process them
        if clade in sranges_map:
            splits_dict = sranges_map[clade]

            for split, support in splits_dict.items():
                # Split has .ancestor and .descendant (both SRangesClade)
                anc_clade = split.ancestor
                desc_clade = split.descendant

                anc_id = clade_to_node_id(anc_clade)
                desc_id = clade_to_node_id(desc_clade)

                # Add child clade nodes (they'll recurse too)
                nodes.append({
                    'data': {
                        'id': anc_id,
                        'label': f"{len(anc_clade.clade) if anc_clade.clade else 0}",
                        'terminal': len(anc_clade.clade) == 1 if anc_clade.clade else False,
                        'range': anc_clade.ancestral_range or ''
                    }
                })

                nodes.append({
                    'data': {
                        'id': desc_id,
                        'label': f"{len(desc_clade.clade) if desc_clade.clade else 0}",
                        'terminal': len(desc_clade.clade) == 1 if desc_clade.clade else False,
                        'range': desc_clade.ancestral_range or ''
                    }
                })

                # Get probabilities
                prob_anc = probability_map.get((clade_id, anc_id), 1.0) if probability_map else 1.0
                prob_desc = probability_map.get((clade_id, desc_id),
                                                1.0) if probability_map else 1.0

                # Add edge to ancestor (avoid duplicates)
                edge_key_anc = (clade_id, anc_id)
                if edge_key_anc not in seen_edges:
                    seen_edges.add(edge_key_anc)
                    edges.append({
                        'data': {
                            'source': clade_id,
                            'target': anc_id,
                            'weight': support,
                            'probability': prob_anc,
                            'role': 'ancestor'
                        }
                    })

                    edges.append({
                        'data': {
                            'source': clade_id,
                            'target': desc_id,
                            'weight': support,
                            'probability': prob_desc,
                            'role': 'descendant'
                        }
                    })

                # Recurse on child clades if they exist as keys
                if anc_clade in sranges_map:
                    add_clade_recursive(anc_clade)
                if desc_clade in sranges_map:
                    add_clade_recursive(desc_clade)

    # Process all clades in the map - handles forests
    for clade in sranges_map.keys():
        if clade not in visited_clades:
            add_clade_recursive(clade)

    # Build HTML template
    html_template = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>SRanges Tree Visualization</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.26.0/cytoscape.min.js"></script>
    <style>
        body {{ margin: 0; padding: 0; overflow: hidden; }}
        #cy {{ width: 100vw; height: 100vh; }}
        #info {{
            position: absolute; top: 10px; left: 10px; z-index: 999;
            background: rgba(0,0,0,0.7); color: white; padding: 10px; 
            border-radius: 5px; max-width: 300px; font-family: monospace;
        }}
        #stats {{
            position: absolute; bottom: 10px; left: 10px; z-index: 999;
            background: rgba(0,0,0,0.7); color: white; padding: 10px;
            border-radius: 5px;
        }}
    </style>
</head>
<body>
    <div id="info">Hover for details</div>
    <div id="stats">Loading...</div>
    <div id="cy"></div>
    <script>
        var nodes = {json.dumps(nodes)};
        var edges = {json.dumps(edges)};

        console.log('Nodes:', nodes.length, 'Edges:', edges.length);
        document.getElementById('stats').innerHTML = 
            '<b>Stats:</b><br>Nodes: ' + nodes.length + 
            '<br>Edges: ' + edges.length;

        var cy = cytoscape({{
            container: document.getElementById('cy'),
            elements: {{ nodes: nodes, edges: edges }},
            style: [
                {{
                    selector: 'node',
                    style: {{
                        'background-color': 'data(terminal) ? "#ff9999" : "#6d4aff"',
                        'border-width': 2,
                        'border-color': '#ffffff',
                        'label': 'data(label)',
                        'font-size': 'data(num_resolutions) ? data(num_resolutions)*2 + 8 : 10px',
                        'color': 'white',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'width': 'data(terminal) ? 25 : 35 + data(num_resolutions, 0)*3',
                        'height': 'data(terminal) ? 25 : 35 + data(num_resolutions, 0)*3'
                    }}
                }},
                {{
                    selector: 'edge',
                    style: {{
                        'width': 'max(1, data(weight))',
                        'line-color': 'rgb(100, ' + Math.floor(data(probability)*150) + ', 100)',
                        'opacity': 'max(0.3, data(probability))',
                        'curve-style': 'bezier',
                        'target-arrow-shape': 'triangle',
                        'arrow-scale': 1.2,
                        'label': 'data(weight)',
                        'font-size': '9px',
                        'color': 'white',
                        'text-outline-color': 'black',
                        'text-outline-width': 1
                    }}
                }},
                {{
                    selector: ':selected',
                    style: {{
                        'border-width': 4,
                        'border-color': '#ffd700',
                        'background-opacity': 1
                    }}
                }},
                {{
                    selector: 'edge:selected',
                    style: {{
                        'line-color': '#ffd700',
                        'width': 4
                    }}
                }}
            ],
            layout: {{ 
                name: 'dagre', 
                rankDir: 'TB',
                nodeSep': 50,
                rankSep': 100
            }},
            minZoom: 0.1,
            maxZoom: 5,
            wheelSensitivity: 0.3,
            selectionType: 'multiselect'
        }});

        // Hover info panel
        cy.on('mouseover', 'node, edge', function(e){{
            var data = e.target.data();
            var html = '<b>' + (data.id || data.role) + '</b><br>';
            if (data.id) {{
                html += 'Terminal: ' + (data.terminal ? 'Yes' : 'No') + '<br>';
                html += 'Range: ' + (data.range || 'None') + '<br>';
                if (data.num_resolutions !== undefined) {{
                    html += 'Resolutions: ' + data.num_resolutions + '<br>';
                }}
            }}
            if (data.weight !== undefined) {{
                html += 'Support: ' + data.weight + '<br>';
                html += 'Probability: ' + data.probability.toFixed(2) + '<br>';
            }}
            document.getElementById('info').innerHTML = html;
        }});

        cy.on('mouseout', 'node, edge', function(e){{
            document.getElementById('info').innerHTML = 'Hover for details<br>Click to select';
        }});

        // Double-click to fit
        cy.on('dblclick', function(){{
            cy.fit();
        }});
    </script>
</body>
</html>'''

    with open(output_file, 'w') as f:
        f.write(html_template)

    print(f"Saved interactive visualization to {output_file}")
    print(f"Graph has {len(nodes)} nodes and {len(edges)} edges")
    return output_file
