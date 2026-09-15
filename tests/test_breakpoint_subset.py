import gzip
import pathlib
import sys
from multiprocessing.pool import ThreadPool

import numpy as np
import pytest
import zstandard as zstd

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "IBDReduce"))
sys.path.insert(0, str(PROJECT_ROOT / "tools"))

import ibdreduce_v3 as reduce
from ibdmap_text_to_parquet import convert


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    # Exercise the CLI pipeline in-process with the same parallel map interface.
    monkeypatch.setattr(reduce.mp, "Pool", ThreadPool)
    gmap = tmp_path / "maps"
    gmap.mkdir()
    with gzip.open(gmap / "test.gmap.gz", "wt") as handle:
        handle.write("pos chrom cM\n")
        for chrom in (1, 2):
            for pos in (1000, 2000, 3000, 4000):
                handle.write(f"{pos} {chrom} {pos / 1000}\n")
    selected = {(1, 1000), (1, 3000), (2, 2000)}
    locations = tmp_path / "locations.txt"
    locations.write_text("# selected sites\nCHROM POS\nchr2 2000\n1 3000\n1 1000\n1 1000\n")
    for chrom in (1, 2):
        for run in (3, 4):
            full, subset = [], []
            for k, pos in enumerate((1000, 2000, 3000, 4000), 1):
                line = f"{chrom} {pos} 0.1 0.2 0.3 {k * chrom / 10} {k / 2**run} {(5-k) / 2**run}\n"
                full.append(line)
                if (chrom, pos) in selected:
                    subset.append(line)
            (tmp_path / f"full_chr{chrom}_run{run}.txt").write_text("".join(full))
            (tmp_path / f"subset_chr{chrom}_run{run}.txt").write_text("".join(subset))

    def run(prefix="full", suffix=".txt", selection=locations, extra=()):
        output = tmp_path / "result.tsv"
        argv = ["ibdreduce_v3.py", "--prefix", str(tmp_path / prefix),
                "--suffix", "_chr{i}_run{j}" + suffix,
                "--nperm", "2", "--nruns", "2", "--at", "3",
                "--gmap", str(gmap), "--output", str(output), *extra]
        if selection is not None:
            argv.extend(["--breakpoints", str(selection)])
        monkeypatch.setattr(sys, "argv", argv)
        reduce.main()
        return output.read_text().splitlines()

    return tmp_path, locations, run


@pytest.mark.parametrize("suffix,extra", [
    (".txt", ()), (".gz", ()), (".zst", ()), (".parquet", ()),
    (".txt", ("--unweighted",)), (".txt", ("--no_avg",)),
    (".txt", ("--two_sided",)), (".txt", ("--fdr",)),
])
def test_subset_matches_physically_filtered_inputs(inputs, suffix, extra):
    root, locations, run = inputs
    for path in root.glob("*_chr*_run*.txt"):
        if suffix == ".gz":
            with gzip.open(path.with_suffix(suffix), "wt") as handle:
                handle.write(path.read_text())
        elif suffix == ".zst":
            with zstd.open(path.with_suffix(suffix), "wt") as handle:
                handle.write(path.read_text())
        elif suffix == ".parquet":
            convert(path, path.with_suffix(suffix), batch_size=2, compression="zstd")

    actual = run(suffix=suffix, extra=extra)
    # Legacy whole-genome operation expects all 22 chromosome files.
    # Empty files on chromosomes 3-22 have no effect on the reference result.
    for chrom in range(3, 23):
        for j in (3, 4):
            (root / f"subset_chr{chrom}_run{j}.txt").write_text("")
    expected = run(prefix="subset", selection=None, extra=extra)
    # Command provenance differs; every result and other metadata must match.
    assert actual[1:] == expected[1:]
    rows = [line.split("\t") for line in actual if not line.startswith("#")][1:]
    assert {(int(row[0]), int(row[1])) for row in rows} == {(1, 1000), (1, 3000), (2, 2000)}
    weights = {(int(row[0]), int(row[1])): float(row[2]) for row in rows}
    assert weights == {(1, 1000): 0.2, (1, 3000): 0.4, (2, 2000): 0.4}


def test_missing_requested_breakpoint_fails(inputs):
    root, locations, run = inputs
    locations.write_text("1 1000\n1 2500\n")
    with pytest.raises(ValueError, match="Requested breakpoints missing"):
        run()


@pytest.mark.parametrize("mutation", ["missing", "reordered"])
def test_inconsistent_later_run_fails(inputs, mutation):
    root, locations, run = inputs
    path = root / "full_chr1_run4.txt"
    lines = path.read_text().splitlines(keepends=True)
    path.write_text("".join(lines[1:] if mutation == "missing" else reversed(lines)))
    with pytest.raises(ValueError, match="does not match|Missing breakpoints"):
        run()


@pytest.mark.parametrize("contents", ["", "# no sites\n", "CHROM POS\n", "1\n", "chrX 100\n", "1 -5\n"])
def test_invalid_breakpoint_list(tmp_path, contents):
    path = tmp_path / "locations.txt"
    path.write_text(contents)
    with pytest.raises(ValueError):
        reduce.read_breakpoints(str(path))


def test_text_skips_unselected_statistics(tmp_path):
    path = tmp_path / "input.txt"
    path.write_text("1 100 invalid statistics\n1 200 0.1 0.2 0.3 0.4 0.5\n")
    rows = list(reduce.iter_result_rows(str(path), {(1, 200)}))
    assert len(rows) == 1
    np.testing.assert_array_equal(rows[0][3], [0.4, 0.5])
