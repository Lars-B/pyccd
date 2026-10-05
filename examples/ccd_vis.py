import json


def sranges_map_to_cytoscape_html(sranges_map, reverse_taxon_map, probability_map=None,
                                  output_file='cytoscape_tree.html'):
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

    def add_clade_recursive(clade, current_depth=0):
        clade_id = clade_to_node_id(clade)

        # Skip if already processed
        if clade_id in visited_clade_ids:
            return

        visited_clade_ids.add(clade_id)
        node_depths[clade_id] = current_depth  # Track depth HERE

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
                prob_anc = probability_map.get((clade_id, anc_id),
                                               1.0) if probability_map else 1.0
                prob_desc = probability_map.get((clade_id, desc_id),
                                                1.0) if probability_map else 1.0

                prob_anc = round(float(prob_anc), 2)
                prob_desc = round(float(prob_desc), 2)

                valid_edges.append({
                    'source': str(clade_id),
                    'target': str(anc_id),
                    'weight': int(support),
                    'probability': prob_anc
                })

                valid_edges.append({
                    'source': str(clade_id),
                    'target': str(desc_id),
                    'weight': int(support),
                    'probability': prob_desc
                })

            # Recurse with incremented depth
            if anc_clade in sranges_map:
                add_clade_recursive(anc_clade, current_depth + 1)
            if desc_clade in sranges_map:
                add_clade_recursive(desc_clade, current_depth + 1)

        # Process all root clades

    for clade in sranges_map.keys():
        add_clade_recursive(clade, 0)

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
            'depth': int(node_depths.get(clade_id, 0))
        })

    # === DEBUG: Verify data integrity ===
    node_ids = {n['id'] for n in nodes}
    orphan_edges = []
    for e in valid_edges:
        if e['source'] not in node_ids or e['target'] not in node_ids:
            orphan_edges.append(e)

    if orphan_edges:
        print(f"WARNING: {len(orphan_edges)} edges reference non-existent nodes!")
        for e in orphan_edges[:3]:
            print(f"  Edge: {e['source']} -> {e['target']}")

    # Remove orphan edges
    valid_edges = [e for e in valid_edges if e['source'] in node_ids and e['target'] in node_ids]

    print(f"Built graph: {len(nodes)} nodes, {len(valid_edges)} edges, {len(all_clades)} clades")

    nodes_json = json.dumps(nodes, separators=(',', ':'))
    edges_json = json.dumps(valid_edges, separators=(',', ':'))

    html_content = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>SRanges Tree</title>
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
    
    data.nodes.sort((a, b) => a.depth - b.depth || b.num_resolutions - a.num_resolutions);
    
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
          'width':2,
          'target-arrow-shape':'triangle',
          'curve-style':'haystack'
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
      alert('ID:'+d.id+'\\nTerminal:'+d.terminal+'\\nLabel:'+d.label);
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
