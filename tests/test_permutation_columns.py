import pathlib
import sys

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.append(str(PROJECT_ROOT / "IBDReduce"))

sys.path.append(str(PROJECT_ROOT / "tools"))

from ibdmap_text_to_parquet import convert
from ibdreduce_v3 import ibdlen, iter_result_rows, parse_line, write_permutation_values


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


def test_parquet_rows_match_result_schema(tmp_path):
    path = tmp_path / "results.parquet"
    table = pa.table({
        "chrom": pa.array(["chr1", "chr1"], type=pa.string()),
        "pos": pa.array([100, 200], type=pa.int32()),
        "orig_cscs_rate": pa.array([0.1, 0.2], type=pa.float64()),
        "orig_cscn_rate": pa.array([0.3, 0.4], type=pa.float64()),
        "orig_cncn_rate": pa.array([0.5, 0.6], type=pa.float64()),
        "original": pa.array([1.25, 1.5], type=pa.float64()),
        "permutations": pa.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]], type=pa.list_(pa.float32())),
    })
    pq.write_table(table, path)

    rows = list(iter_result_rows(str(path)))

    assert rows[0][0:2] == (1, 100)
    np.testing.assert_allclose(rows[0][2], np.array([0.1, 0.3, 0.5]))
    np.testing.assert_allclose(rows[0][3], np.array([1.25, 0.1, 0.2, 0.3]), rtol=1e-6)
    assert rows[1][0:2] == (1, 200)
    np.testing.assert_allclose(rows[1][2], np.array([0.2, 0.4, 0.6]))
    np.testing.assert_allclose(rows[1][3], np.array([1.5, 0.4, 0.5, 0.6]), rtol=1e-6)


def test_parquet_ibdlen_reads_positions(tmp_path):
    class DummyMap:
        gmap = {1: {100: 1.0, 200: 2.5}}

    path = tmp_path / "results.parquet"
    table = pa.table({
        "chrom": pa.array(["chr1", "chr1"], type=pa.string()),
        "pos": pa.array([100, 200], type=pa.int32()),
        "orig_cscs_rate": pa.array([0.1, 0.2], type=pa.float64()),
        "orig_cscn_rate": pa.array([0.3, 0.4], type=pa.float64()),
        "orig_cncn_rate": pa.array([0.5, 0.6], type=pa.float64()),
        "original": pa.array([1.25, 1.5], type=pa.float64()),
        "permutations": pa.array([[0.1, 0.2], [0.3, 0.4]], type=pa.list_(pa.float32())),
    })
    pq.write_table(table, path)

    breakpoints, total = ibdlen(str(path), DummyMap())

    assert breakpoints == 2
    assert total == 2.5


def test_convert_text_output_to_parquet(tmp_path):
    text_path = PROJECT_ROOT / "tests" / "data" / "sample_chr1_run0.txt"
    parquet_path = tmp_path / "converted.parquet"

    rows_written = convert(text_path, parquet_path, batch_size=2, compression="zstd")

    assert rows_written == 3
    rows = list(iter_result_rows(str(parquet_path)))
    assert rows[0][0:2] == (1, 1000)
    np.testing.assert_allclose(rows[0][2], np.array([0.1, 0.2, 0.3]))
    np.testing.assert_allclose(rows[0][3], np.array([0.3, 0.4, 0.1]), rtol=1e-6)
    assert rows[-1][0:2] == (1, 3000)
