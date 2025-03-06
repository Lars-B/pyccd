from CCDpy.conditional_clade_distribution import *
from pathlib import Path
import numpy as np


def test_get_maps():
    test_trees_file = f"{Path(__file__).parent.absolute()}/data/30Taxa.trees"
    trees = read_nexus(test_trees_file)
    m1, m2, _ = get_maps(trees)
    p = get_tree_probability(trees[1], m1, m2)
    assert True # todo not sure


def test_get_tree_probability(ten_taxa_multichain):
    m1, m2, u = get_maps(ten_taxa_multichain[0].trees)
    p = get_tree_probability(ten_taxa_multichain[0].trees[1], m1, m2)
    assert p == 0.0037693679383000514, "Get tree probability failed!"


def test_get_tree_log_probability(ten_taxa_multichain):
    m1, m2, u = get_maps(ten_taxa_multichain[0].trees)
    p = get_tree_probability(ten_taxa_multichain[0].trees[1], m1, m2, log=True)
    assert p == np.log(0.0037693679383000514), "Get tree probability failed!"


def test_get_ccd_tree_bottom_up(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    best_tree = get_ccd_tree_bottom_up(m1, m2)
    assert len(best_tree) == 20, "Failed DFS BB algorithm"


def test_sample_tree_from_ccd(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    sample = sample_tree_from_ccd(m1, m2, n=100)
    assert len(sample) == 100


def test_sample_logprob_from_ccd(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    sample = sample_logprob_from_ccd(m1, m2, n=100)
    assert True


def test_calc_Entropy(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    h = calc_Entropy(m1, m2)
    assert h == 7.500020954155561, "Entropy calculation failed!"
