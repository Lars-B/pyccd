import copy
from collections import Counter

import numpy as np
from tqdm import tqdm

from brokilon.ccd.domain.sranges import sranges


def smallest_common_superclade(s, tree_clades, rare):
    """Smallest non-rare clade in the same tree whose taxa include all of s's taxa."""
    candidates = [
        c for c in tree_clades
        if c not in rare and c != s and s.clade <= c.clade
    ]
    return min(candidates, key=lambda c: len(c.clade), default=None)


def flag_outliers(scores, z_cutoff=1.0):
    """Indices of trees whose log score is more than z_cutoff robust deviations above
    the median, with the median and MAD taken from trees scoring above 0."""
    scores = np.asarray(scores, dtype=float)
    log_scores = np.log1p(scores)
    positive = log_scores[scores > 0]
    print(f"{len(positive)} of {len(scores)} trees have a score above 0")

    if len(positive) == 0:
        print("No tree has a score above 0, no trees flagged.")
        return np.array([], dtype=int)

    med = np.median(positive)
    mad = np.median(np.abs(positive - med)) * 1.4826
    if mad == 0:
        # Most scored trees share one value, so flag only the ones above it
        flagged = np.where(log_scores > med)[0]
    else:
        robust_z = (log_scores - med) / mad
        flagged = np.where(robust_z > z_cutoff)[0]

    print(f"{len(flagged)} of {len(scores)} trees flagged (z > {z_cutoff})")
    return flagged


def prune_outlier_trees(trees, taxon_map, k=1, z_cutoff=1.0):
    # Clades of each tree, each from its own fresh copy
    per_tree_clades = []
    for tree in tqdm(trees, desc="Collecting clades", unit="tree"):
        counts, _, _, _ = sranges.get_sranges_map([copy.deepcopy(tree)], taxon_map)
        per_tree_clades.append(set(counts))

    # Count across trees by combining the per-tree results
    clade_counts = Counter()
    for clades in per_tree_clades:
        clade_counts.update(clades)

    rare = {c for c, n in clade_counts.items() if n <= k}

    # For each rare clade: (rare clade size, superclade size,
    #                       rare clade count, superclade count)
    rows = []
    for tc in per_tree_clades:
        tree_rows = []
        for s in tc:
            if s not in rare:
                continue
            sup = smallest_common_superclade(s, tc, rare)
            tree_rows.append((
                len(s.clade),
                len(sup.clade) if sup else None,
                clade_counts[s],
                clade_counts[sup] if sup else None,
            ))
        rows.append(tree_rows)

    scores = []
    for tree_rows in rows:
        size_score = 0.0
        count_score = 0.0
        weighted_count_score = 0.0
        for len_s, len_sup, count_s, count_sup in tree_rows:
            if len_sup is None:
                continue
            size_score += len_sup / len_s
            count_score += count_sup / count_s
            weighted_count_score += (len_sup - len_s) * count_sup / count_s
        scores.append((size_score, count_score, weighted_count_score))

    # for name, j in [("size", 0), ("count", 1), ("weighted count", 2)]:
    #     ranked = sorted(range(len(scores)), key=lambda i: scores[i][j], reverse=True)
    #     print(f"Top 10 trees by {name} score:")
    #     for i in ranked[:10]:
    #         print(f"  tree {i}: {scores[i][j]:.2f}  ({len(rows[i])} rare clades)")

    weighted = [s[2] for s in scores]
    flagged = set(flag_outliers(weighted, z_cutoff=z_cutoff).tolist())

    kept = [
        copy.deepcopy(tree)
        for i, tree in enumerate(trees)
        if i not in flagged
    ]
    # print(f"Kept {len(kept)} of {len(trees)} trees.")
    return kept
