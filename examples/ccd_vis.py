import json


def sranges_map_to_cytoscape_html(sranges_map, reverse_taxon_map, clade_count_map,
                                  output_file='sranges-ccd-vis.html'):
    """
    Minimal, validated Cytoscape HTML generator.
    sranges_map: Dict[SRangesClade, Dict[AncestralSplit, int]]
    """

    # === PHASE 1: Collect unique clades ===
    all_clades = {}  # clade_id → SRangesClade

    def clade_to_node_id(clade):
        """Generate stable ID even for empty frozensets."""
        nonlocal reverse_taxon_map
        # Case 1: Has taxa - use sorted taxa + range
        range_id = str(
            reverse_taxon_map[clade.ancestral_range]) if clade.ancestral_range else 'NR'
        if clade.clade:
            taxa = tuple(sorted(clade.clade))
            # base_hash = hash(taxa)
            clade_string = str({reverse_taxon_map[t] for t in taxa})
        # Case 2: Empty frozenset (special/range node) - use object id + range
        else:
            # Use object id as fallback for uniqueness
            obj_id = id(clade)
            # base_hash = obj_id
            clade_string = obj_id

        # Combine into stable ID
        clean_id = f"c_{clade_string}_{range_id}".replace(
            '"', '').replace('\n', '').strip()

        # Final safety: ensure non-empty ID
        # if not clean_id or len(clean_id) < 5:
        #     clean_id = f"clade_special_{}"

        return str(clean_id)

    # === PHASE 2: Traverse, collect edges AND track depths ===
    visited_clade_ids = set()
    valid_edges = []
    node_depths = {}  # clade_id → depth integer
    seen_edge_pairs = set()

    def add_clade_recursive(clade):

        clade_id = clade_to_node_id(clade)

        # Skip if already processed
        if clade_id in visited_clade_ids:
            return

        visited_clade_ids.add(clade_id)
        node_depths[clade_id] = len(clade)  # Track depth HERE

        # Register this clade
        all_clades[clade_id] = clade

        # Get splits for this clade
        splits_dict = sranges_map.get(clade, {})

        for split, support in splits_dict.items():
            anc_clade = split.ancestor
            desc_clade = split.descendant

            anc_id = clade_to_node_id(anc_clade)
            desc_id = clade_to_node_id(desc_clade)

            # Validate
            if anc_id not in all_clades:
                all_clades[anc_id] = anc_clade
            if desc_id not in all_clades:
                all_clades[desc_id] = desc_clade

            # Create edges - VALIDATE SOURCE/TARGET
            if anc_id and desc_id and clade_id:
                e1 = (clade_id, anc_id)
                e2 = (clade_id, desc_id)

                if e1 not in seen_edge_pairs:
                    seen_edge_pairs.add(e1)
                    valid_edges.append({
                        'source': str(clade_id),
                        'target': str(anc_id),
                        'weight': int(support),
                    })
                if e2 not in seen_edge_pairs:
                    seen_edge_pairs.add(e2)
                    valid_edges.append({
                        'source': str(clade_id),
                        'target': str(desc_id),
                        'weight': int(support),
                    })

            # Recurse with incremented depth
            if anc_clade in sranges_map:
                add_clade_recursive(anc_clade)
            if desc_clade in sranges_map:
                add_clade_recursive(desc_clade)

    # Process all root clades
    for clade in sranges_map.keys():
        add_clade_recursive(clade)

    # === Now build edges (existing code) ===

    # === PHASE 3: Build final node/edge arrays ===
    nodes = []
    for clade_id, clade in all_clades.items():
        is_terminal = len(clade.clade) == 1 if clade.clade else False
        num_resolutions = len(sranges_map.get(clade, {}))

        nodes.append({
            'id': str(clade_id),
            'label': str(len(clade.clade)) if clade.clade else '0',
            'terminal': bool(is_terminal),
            'num_resolutions': int(num_resolutions),
            'freq': clade_count_map.get(clade, -1),
        })

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
#debug{{position:absolute;top:10px;left:10px;background:#333;color:#fff;padding:10px;font-family:monospace;z-index:999;max-width:300px;}}
</style>
</head>
<body>
<div id="debug">Loading...</div>
<div id="cy"></div>
<script>
window.onload = function() {{
  try {{
    const data = {{nodes:{nodes_json},edges:{edges_json}}};
    const elements = {{
      nodes: data.nodes.map(n => ({{ data: n }})),
      edges: data.edges.map(e => ({{ data: e }}))
    }};
    // Debug output
    document.getElementById('debug').innerHTML = 
      'Nodes: ' + data.nodes.length + '<br>' +
      'Edges: ' + data.edges.length + '<br>' +
      'Status: Building...';
    
    const cy = cytoscape({{
      container: document.getElementById('cy'),
      elements: elements,
      style: [
        {{selector:'node',style:{{
          'background-color':'#6d4aff',
          'border-color':'#ffffff',
          'border-width':2,
          'label':'data(id)',
          'color':'#000000',
          'font-size':10,
          'text-valign':'center',
          'text-halign':'center',
          'width':'mapData(num_resolutions, 0, 20, 35, 80)',
          'height':'mapData(num_resolutions, 0, 20, 35, 80)',
        }}}},
        {{selector:'node[terminal=true]',style:{{
          'background-color':'#ff9999'
        }}}},
        {{selector:'edge',style:{{
          'line-color':'#64bf64',
          'width':4,
          'curve-style': 'bezier',
          'target-arrow-shape':'triangle',
          'arrow-scale': 2.0,
          'targer-arrow-color': "#ccc"
        }}}},
        {{selector:':selected',style:{{
          'border-color':'#ffd700',
          'border-width':4
        }}}}
      ],
      layout: {{
            name: 'breadthfirst',
            directed: true,
            padding: 100,
            spacingFactor: 2.5,  // More spread between levels
            avoidOverlap: 0.8,   // Prevent node overlap
            pack: true,          // Handle disconnected components
            roots: function(node) {{
                    // Put high-depth nodes at top (invert order)
                    return node.data('depth') === 0;
                    }}
        }},
      minZoom:0.1,
      maxZoom:3
    }});

    document.getElementById('debug').innerHTML = 
      'Nodes: ' + data.nodes.length + '<br>' +
      'Edges: ' + data.edges.length + '<br>' +
      'Status: Ready - hover/click nodes';

    // Click handler
    cy.on('tap','node',function(evt){{
      const nd=evt.target;
      const d=nd.data();
      alert('ID:'+d.id+'\\nTerminal:'+d.terminal+'\\nLabel:'+d.label+'\\nFreq:'+d.freq);
    }});

    console.log('Cytoscape initialized successfully');
  }} catch(e) {{
    document.getElementById('debug').innerHTML = 'ERROR: '+e.message;
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
