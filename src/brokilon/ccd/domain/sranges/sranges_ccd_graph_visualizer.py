import json

FONT_SIZE = 30


def label_width(text, char_px=8, padding=2, min_w=30):
    """Estimate the minimum width required to display a node label inside the node.

    Args:
        text: The label text.
        char_px: Approximate pixel width per character (currently unused).
        padding: Additional horizontal padding in pixels.
        min_w: Minimum node width in pixels.

    Returns:
        The estimated label width in pixels.
    """
    return max(min_w, int(len(text) * FONT_SIZE * 0.8 + padding))


def sranges_map_to_cytoscape_html(sranges_map, reverse_taxon_map, clade_count_map,
                                  output_file='sranges-ccd-vis.html', use_reverse_taxon_map=True):
    """Generate an interactive Cytoscape HTML visualization of the stratigraphic-range CCD graph.

    Args:
        sranges_map: Mapping from stratigraphic-range clades (SRangesClade) to
            mappings of ancestral splits (AncestralSplit) and their counts.
        reverse_taxon_map: Mapping from internal taxon identifiers to their
            original taxon labels.
        clade_count_map: Mapping from stratigraphic-range clades to their counts.
        output_file: Path where the generated HTML visualization is saved.
            Defaults to 'sranges-ccd-vis.html'.
        use_reverse_taxon_map: Whether to use the reverse taxon map to display
            original taxon labels instead of internal identifiers. Defaults to True.

    Returns:
        None. Writes the interactive Cytoscape visualization to an HTML file.
    """

    # === PHASE 1: Collect unique clades ===
    all_clades = {}  # clade_id → SRangesClade

    def clade_to_node_id(clade):
        """Generate stable ID even for empty frozensets."""

        def taxon_label(t):
            if use_reverse_taxon_map:
                return str(reverse_taxon_map[t])
            return str(t).replace("_first", "")

        # Case 1: Has taxa - use sorted taxa + range
        range_id = taxon_label(clade.ancestral_range) if clade.ancestral_range else '∅'
        if clade.clade:
            taxa = sorted(clade.clade)
            clade_string = '{' + ', '.join(sorted(taxon_label(t) for t in taxa)) + '}'
        # Case 2: Empty frozenset (special/range node) - use object id + range
        else:
            clade_string = '∅'

        # Combine into stable ID
        clean_id = f"{clade_string}, {range_id}".replace(
            '"', '').replace('\n', '').strip()

        return str(clean_id)

    # === PHASE 2: Traverse, collect edges AND track depths ===
    visited_clade_ids = set()
    valid_edges = []
    split_nodes = []
    seen_split_ids = set()
    child_edge_map = {}

    def add_clade_recursive(clade):

        nonlocal sranges_map, clade_count_map

        clade_id = clade_to_node_id(clade)

        # Skip if already processed
        if clade_id in visited_clade_ids:
            return

        visited_clade_ids.add(clade_id)

        # Register this clade
        all_clades[clade_id] = clade

        # Get splits for this clade
        splits_dict = sranges_map.get(clade, {})

        for split, split_count in splits_dict.items():
            anc_clade = split.ancestor
            desc_clade = split.descendant

            anc_id = clade_to_node_id(anc_clade)
            desc_id = clade_to_node_id(desc_clade)

            if anc_id not in all_clades:
                all_clades[anc_id] = anc_clade
            if desc_id not in all_clades:
                all_clades[desc_id] = desc_clade

            # One shared split node per unique (ancestor, descendant) pair
            split_id = f"split|{anc_id}|{desc_id}"
            split_label = f"{anc_id} || {desc_id}"

            if split_id not in seen_split_ids:
                seen_split_ids.add(split_id)
                split_nodes.append({
                    'id': split_id,
                    'label': split_label,
                    'is_split': True,
                    'width': label_width(split_label)
                })

            # Parent clade -> split: one edge per parent, annotation is per-parent
            valid_edges.append({
                'source': str(clade_id),
                'target': split_id,
                'to_split': True,
                'ccd_prob': f"{split_count} / {clade_count_map[clade]}"
            })

            # Split -> child clades: one edge per (split, child), counts summed
            for target in {str(anc_id), str(desc_id)}:
                edge = child_edge_map.setdefault(
                    (split_id, target),
                    {'source': split_id, 'target': target, 'count': 0}
                )
                edge['count'] += split_count

            if anc_clade in sranges_map:
                add_clade_recursive(anc_clade)
            if desc_clade in sranges_map:
                add_clade_recursive(desc_clade)

    # Process all root clades
    for clade in sranges_map.keys():
        add_clade_recursive(clade)

    # Add the merged split -> child edges once everything has been counted
    valid_edges.extend(child_edge_map.values())

    # === Now build edges (existing code) ===

    # === PHASE 3: Build final node/edge arrays ===
    nodes = []
    for clade_id, clade in all_clades.items():
        is_terminal = len(clade.clade) == 1 if clade.clade else False
        num_resolutions = len(sranges_map.get(clade, {}))
        clade_label = str(len(clade.clade)) if clade.clade else '0'
        nodes.append({
            'id': str(clade_id),
            'label': clade_label,
            'width': label_width(clade_id),
            'num_resolutions': int(num_resolutions),
            'freq': clade_count_map.get(clade, -1),
        })
    nodes.extend(split_nodes)

    print(f"Built graph: {len(nodes)} nodes, {len(valid_edges)} edges, {len(all_clades)} clades")

    nodes_json = json.dumps(nodes, separators=(',', ':'))
    edges_json = json.dumps(valid_edges, separators=(',', ':'))

    html_content = f'''<!DOCTYPE html>
    <html>
    <head>
    <meta charset="UTF-8">
    <title>SRanges CCD</title>
    <script src="https://cdn.jsdelivr.net/npm/cytoscape@3.26.0/dist/cytoscape.min.js"></script>
    <style>
    #cy{{width:100vw;height:100vh;background-color:#e5e5e5;}}
    #controls{{position:absolute;top:10px;left:10px;background:#333;color:#fff;padding:10px;border-radius:5px;font-family:monospace;z-index:999;max-width:300px;}}
    #search{{width:100%;padding:5px;margin-bottom:5px;font-size:12px;box-sizing:border-box;}}
    #results{{margin-top:5px;font-size:11px;color:#aaa;}}
    #info{{position:absolute;top:70px;left:10px;background:#333;color:#fff;padding:10px;border-radius:5px;font-family:monospace;max-width:300px;z-index:999;}}
    #stats{{position:absolute;bottom:10px;left:10px;background:#333;color:#fff;padding:10px;border-radius:5px;font-family:monospace;z-index:999;}}
    </style>
    </head>
    <body>
    <div id="controls">
        <input type="text" id="search" placeholder="Search clade by ID, label, or range...">
        <div id="results">0 matches</div>
    </div>
    <div id="info">Hover or click nodes for details</div>
    <div id="stats">Loading graph...</div>
    <div id="cy"></div>
    <script>
    window.onload = function() {{
      try {{
        const data = {{nodes:{nodes_json},edges:{edges_json}}};
        const elements = {{
          nodes: data.nodes.map(n => ({{ data: n }})),
          edges: data.edges.map(e => ({{ data: e }}))
        }};

        document.getElementById('stats').innerHTML = 
          'Nodes: ' + data.nodes.length + '<br>' +
          'Edges: ' + data.edges.length + '<br>' +
          'Status: Building...';

        const cy = cytoscape({{
          container: document.getElementById('cy'),
          elements: elements,
          style: [
            {{selector:'node',style:{{
            'shape': 'round-rectangle',
              'background-color':'#b8c7ff',
              'border-color':'#ffffff',
              'border-width':2,
              'label':'data(id)',
              'color':'#000000',
              'font-size':{FONT_SIZE},
              'text-valign':'center',
                'text-halign':'center',
                'text-wrap':'wrap',
                'text-max-width':'data(width)',
                'width': 'data(width)',
                'height': '30',
            }}}},
            {{selector:'node[?is_split]',style:{{
                'shape':'rectangle',
                'background-color':'#a8e6a1',
                'label':'data(label)',
            }}}},
            {{selector:'edge',style:{{
                'line-color':'#64bf64',
                'width':4,
                'curve-style':'bezier',
                'target-arrow-shape':'triangle',
                'arrow-scale':2.0,
                'target-arrow-color':'#64bf64',
                'label':'data(count)',
                'font-size':{FONT_SIZE},
                'color':'#000000'
            }}}},
            {{selector:'edge[?to_split]',style:{{
                'label':'data(ccd_prob)',
            }}}},
            {{selector:':selected',style:{{
              'border-color':'#ffd700',
              'border-width':4
            }}}}
          ],
          layout: {{
                name: 'breadthfirst',
                fit: 'true',
                directed: true,
                padding: 10,
                spacingFactor: 0.8,
                avoidOverlap: 1.2,
                pack: true,
                roots: function(node) {{
                  return node.data('depth') === 0;
                }}
            }},
          minZoom:0.1,
          maxZoom:3
        }});

        // Initialize search
        const searchInput = document.getElementById('search');
        const resultsDiv = document.getElementById('results');

        searchInput.addEventListener('keyup', function() {{
          const query = this.value.toLowerCase().trim();
          const matches = [];

          cy.nodes().removeClass('matched');

          if (query === '') {{
            resultsDiv.innerText = '0 matches';
            return;
          }}

          cy.nodes().forEach(node => {{
            const d = node.data();
            const searchText = [
              d.id || '',
              d.label || '',
              d.range || '',
              String(d.terminal || '')
            ].join(' ').toLowerCase();

            if (searchText.includes(query)) {{
              node.addClass('matched');
              matches.push(node);
            }}
          }});

          if (matches.length > 0) {{
            const firstMatch = cy.nodes('.matched').first();
            firstMatch.select();
            cy.fit(firstMatch, 100);
          }}

          resultsDiv.innerText = matches.length + ' match' + (matches.length !== 1 ? 'es' : '') + ' found';
        }});

        // Click handler
        cy.on('tap','node',function(evt){{
          const nd=evt.target;
          const d=nd.data();
          const info = 
            '<b>ID:</b> ' + d.id + '<br>' +
            '<b>Terminal:</b> ' + d.terminal + '<br>' +
            '<b>Label:</b> ' + d.label + '<br>' +
            '<b>Depth:</b> ' + (d.depth || 'N/A') + '<br>' +
            '<b>Resolutions:</b> ' + (d.num_resolutions || 'N/A');
          document.getElementById('info').innerHTML = info;
        }});

        // Hover handler
        cy.on('mouseover','node',function(evt){{
          const d=evt.target.data();
          document.getElementById('info').innerHTML = 
            '<b>ID:</b> ' + d.id + '<br>' +
            '<b>Terminal:</b> ' + d.terminal + '<br>' +
            '<b>Label:</b> ' + d.label + '<br>' +
            '<b>Depth:</b> ' + (d.depth || 'N/A') + '<br>' +
            '<b>Resolutions:</b> ' + (d.num_resolutions || 'N/A');
        }});

        cy.on('mouseout','node',function(){{
          document.getElementById('info').innerHTML = 'Hover or click nodes for details';
        }});

        document.getElementById('stats').innerHTML = 
          'Nodes: ' + data.nodes.length + '<br>' +
          'Edges: ' + data.edges.length + '<br>' +
          'Status: Ready - use search box above';

        console.log('Cytoscape initialized successfully');
      }} catch(e) {{
        document.getElementById('stats').innerHTML = 'ERROR: '+e.message;
        console.error(e);
      }}
    }};
    </script>
    </body>
    </html>'''

    with open(output_file, 'w') as f:
        f.write(html_content)

    print(f"Saved: {output_file}")
    return output_file
