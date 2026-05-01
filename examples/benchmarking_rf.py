import time
from pathlib import Path

from brokilon.core import read_nexus_trees
from brokilon.metrics import robinson_foulds_ccd, robinson_foulds

tree_file = f"{Path(__file__).parent.absolute().parent}/examples/data/400Taxa.trees"
trees = read_nexus_trees(tree_file)

t1 = trees[10]
t2 = trees[20]


def benchmark(func, n_repeats=10):
    start = time.perf_counter()
    for _ in range(n_repeats):
        func(t1, t2)
    end = time.perf_counter()
    return (end - start) / n_repeats


time_v1 = benchmark(robinson_foulds, n_repeats=10)
time_v2 = benchmark(robinson_foulds_ccd, n_repeats=10)

print(f"Average time for regular RF: {time_v1}")
print(f"Average time for ccd RF: {time_v2}")
print(f"Ratio of timings is {time_v1 / time_v2}")
