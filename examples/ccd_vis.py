import json


def sranges_map_to_cytoscape_html(sranges_map, probability_map=None,
                                  output_file='cytoscape_tree.html'):
    """
    Minimal, validated Cytoscape HTML generator.
    sranges_map: Dict[SRangesClade, Dict[AncestralSplit, int]]
    """

    # === PHASE 1: Collect unique clades ===
    all_clades = {}  # clade_id → SRangesClade

    def clade_to_node_id(clade):
        """Generate stable ID even for empty frozensets."""

        # Case 1: Has taxa - use sorted taxa + range
        if clade.clade:
            taxa = tuple(sorted(clade.clade))
            range_id = str(clade.ancestral_range) if clade.ancestral_range else 'NO_RANGE'
            base_hash = hash(taxa)
        # Case 2: Empty frozenset (special/range node) - use object id + range
        else:
            # Use object id as fallback for uniqueness
            obj_id = id(clade)
            range_id = str(clade.ancestral_range) if clade.ancestral_range else 'NO_RANGE'
            base_hash = obj_id

        # Combine into stable ID
        clean_id = f"clade_{base_hash % 10000}_{range_id}".replace('"', '').replace('\n',
                                                                                    '').strip()

        # Final safety: ensure non-empty ID
        if not clean_id or len(clean_id) < 5:
            clean_id = f"clade_special_{base_hash}"

        return str(clean_id)

    # === PHASE 2: Traverse and collect edges ===
    visited_clade_ids = set()
    valid_edges = []  # List of dicts with source, target, weight

    def add_clade_recursive(clade):
        clade_id = clade_to_node_id(clade)

        # Skip if already processed
        if clade_id in visited_clade_ids:
            return

        visited_clade_ids.add(clade_id)

        # Register this clade
        all_clades[clade_id] = clade

        # Get splits for this clade
        splits_dict = sranges_map.get(clade, {})

        for split, support in splits_dict.items():
            anc_clade = split.ancestor
            desc_clade = split.descendant

            anc_id = clade_to_node_id(anc_clade)
            desc_id = clade_to_node_id(desc_clade)

            # Validate: both child clades must exist as keys in sranges_map OR be terminal
            # Either way, register them so they appear as nodes
            if anc_id not in all_clades:
                all_clades[anc_id] = anc_clade
            if desc_id not in all_clades:
                all_clades[desc_id] = desc_clade

            # Create edges - VALIDATE SOURCE/TARGET
            if anc_id and desc_id and clade_id:
                prob_anc = probability_map.get((clade_id, anc_id), 1.0) if probability_map else 1.0
                prob_desc = probability_map.get((clade_id, desc_id),
                                                1.0) if probability_map else 1.0

                # Round probabilities to avoid floating point weirdness
                prob_anc = round(float(prob_anc), 2)
                prob_desc = round(float(prob_desc), 2)

                # if not anc_id or not desc_id or not clade_id:
                #     print(f"SKIPPING EDGE: parent={clade_id}, anc={anc_id}, desc={desc_id}")
                #     continue  # Skip this split

                print(str(clade_id))

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

            # Recurse on children if they have splits
            if anc_clade in sranges_map:
                add_clade_recursive(anc_clade)
            if desc_clade in sranges_map:
                add_clade_recursive(desc_clade)

    # Process all root clades
    for clade in sranges_map.keys():
        add_clade_recursive(clade)

    # === PHASE 3: Build final node/edge arrays ===
    nodes = []
    for clade_id, clade in all_clades.items():
        is_terminal = len(clade.clade) == 1 if clade.clade else False
        num_resolutions = len(sranges_map.get(clade, {}))

        nodes.append({
            'id': str(clade_id),
            'label': str(len(clade.clade)) if clade.clade else '0',
            'terminal': bool(is_terminal),
            'num_resolutions': int(num_resolutions)
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

    # === PHASE 4: Generate minimal HTML ===

    # print("\n=== EDGE VALIDATION ===")
    # for i, e in enumerate(valid_edges):
    #     src = e.get('source')
    #     tgt = e.get('target')
    #     if not src or src == 'None' or src == 'null' or src == '':
    #         print(f"BAD EDGE {i}: source={repr(src)}, target={repr(tgt)}")
    #     if not tgt or tgt == 'None' or tgt == 'null' or tgt == '':
    #         print(f"BAD EDGE {i}: source={repr(src)}, target={repr(tgt)}")
    # print("=" * 25)

    nodes_json = json.dumps(nodes, separators=(',', ':'))
    edges_json = json.dumps(valid_edges, separators=(',', ':'))



    html_content = f'''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>SRanges Tree</title>
<script src="https://cdn.jsdelivr.net/npm/cytoscape@3.26.0/dist/cytoscape.min.js"></script>
<style>
#cy{{width:100vw;height:100vh}}
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
          'label':'data(label)',
          'color':'#ffffff',
          'font-size':10,
          'text-valign':'center',
          'text-halign':'center',
          'width':35,
          'height':35
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
      layout:{{name:'grid',rows:Math.ceil(Math.sqrt(data.nodes.length))}},
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
