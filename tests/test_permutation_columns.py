import pathlib
import sys

import numpy as np


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.append(str(PROJECT_ROOT / "IBDReduce"))

from ibdreduce_v3 import parse_line, write_permutation_values


def test_permutation_columns_do_not_overlap():
    data_dir = PROJECT_ROOT / "tests" / "data"
    run_paths = [data_dir / f"sample_chr1_run{j}.txt" for j in range(2)]

    first_line = run_paths[0].read_text(encoding="utf-8").splitlines()[0]
    _, _, _, first_vals = parse_line(first_line)
    nperm = len(first_vals) - 1

    with run_paths[0].open("r", encoding="utf-8") as handle:
        breakpoints = sum(1 for _ in handle)

    data = np.zeros((breakpoints, nperm * len(run_paths) + 1), dtype=float)

    for offset, path in enumerate(run_paths):
        idx = 0
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                _, _, _, vals = parse_line(line)
                write_permutation_values(data, idx, offset, vals, nperm)
                idx += 1

    for offset, path in enumerate(run_paths):
        start = 1 + offset * nperm
        end = 1 + (offset + 1) * nperm
        idx = 0
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                _, _, _, vals = parse_line(line)
                np.testing.assert_allclose(data[idx, start:end], vals[1:])
                idx += 1
