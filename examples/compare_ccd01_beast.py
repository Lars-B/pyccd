import subprocess
from pathlib import Path

from pyccd import read_nexus_trees
from pyccd.ccd0_attempt import get_ccd0, get_tree_probability, get_map_tree
from pyccd.ccd import get_tree_probability as get_ccd1_tree_probability, get_ccd_tree_bottom_up
from pyccd.ccd import get_maps
from pyccd.tree import TreeNode


def compare_to_java_tree_probs(java_tree_probs, python_tree_probs):
    tolerance = 1e-12

    prob_problem = []
    for k in python_tree_probs.keys():
        v1 = python_tree_probs[k]
        v2 = java_tree_probs[k]
        if abs(v1 - v2) > tolerance:
            prob_problem.append((k, v1, v2))
    if len(prob_problem) > 0:
        print(prob_problem)
    else:
        print("No mismatches outside of tolerance were found!")


def parse_jave_output(java_output):
    java_tree_probs = {}
    csv_lines = java_output.split("treeIndex,probability\n", 1)[1].splitlines()
    for line in csv_lines:
        ix, prob = line.split(",")
        java_tree_probs[int(ix)] = float(prob)
    return java_tree_probs


def compare_tree_probs():
    tree_files = [
        f"{Path(__file__).parent.absolute().parent}/examples/data/30Taxa.trees",
        # f"{Path(__file__).parent.absolute().parent}/examples/data/roetzer40.trees"
    ]

    for ccd_type in ["CCD0", "CCD1"]:
        for tree_file in tree_files:
            cmd = [
                "/Applications/BEAST 2.7.7/bin/applauncher",
                "CCDProbCalculator",
                "-trees", tree_file,
                "-burnin", "0",
                "-ccdType", ccd_type
            ]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            java_tree_probs = parse_jave_output(result.stdout)

            python_tree_probs = {}
            trees = read_nexus_trees(tree_file, breath_trees=False, label_transm_history=False)
            match ccd_type:
                case "CCD0":
                    partitions_ccp = get_ccd0(trees)

                    for i, t in enumerate(trees):
                        python_tree_probs[i] = get_tree_probability(t, partitions_ccp)
                case "CCD1":
                    m1, m2, _ = get_maps(trees)
                    for i, t in enumerate(trees):
                        python_tree_probs[i] = get_ccd1_tree_probability(t, m1, m2)
                case _:
                    raise ValueError("Unsupported...")

            print(f"Evaluating {tree_file} with type {ccd_type}")
            compare_to_java_tree_probs(java_tree_probs, python_tree_probs)


def compare_map_trees():
    tree_files = [
        f"{Path(__file__).parent.absolute().parent}/examples/data/30Taxa.trees",
        # f"{Path(__file__).parent.absolute().parent}/examples/data/roetzer40.trees"
    ]

    for ccd_type in ["CCD0", "CCD1"]:
        for tree_file in tree_files:
            cmd = [
                "/Applications/BEAST 2.7.7/bin/treeannotator",
                "-burnin", "0",
                "-topology", ccd_type,
                "-file", tree_file,
            ]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            nexus_string = "#NEXUS" + result.stdout.split("#NEXUS")[1]

            import tempfile

            with tempfile.NamedTemporaryFile(mode="w", suffix=".trees", delete=True) as f:
                f.write(nexus_string)
                f.flush()
                java_tree = read_nexus_trees(f.name, breath_trees=False,
                                         label_transm_history=False)

            trees = read_nexus_trees(tree_file, breath_trees=False, label_transm_history=False)
            match ccd_type:
                case "CCD0":
                    partitions_ccp = get_ccd0(trees)
                    map_tree = get_map_tree(partitions_ccp)
                case "CCD1":
                    m1, m2, _ = get_maps(trees)
                    map_tree_str = get_ccd_tree_bottom_up(m1, m2)
                    map_tree = TreeNode(newick=map_tree_str)
                case _:
                    raise ValueError("Unsupported...")
            dist = map_tree.robinson_foulds(java_tree[0])
            if dist[0] != 0:
                print(f"ERROR: Map tree difference detected for {ccd_type}: {dist[0]} "
                      f"({tree_file})")
            else:
                print(f"No difference in MAP tree topology for {ccd_type}: {dist[0]} "
                      f"({tree_file})")


if __name__ == '__main__':
    # compare_tree_probs()
    compare_map_trees()
