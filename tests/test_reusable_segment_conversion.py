from collections import Counter
import gzip
import mmap
from pathlib import Path
import subprocess
import sys

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "tools"))

from compare_breakpoint_loading import decode_block, read_header, read_index
from ibd_to_binary import BLOCK_CHECKPOINT, Indices, read_segments


@pytest.mark.parametrize("script", ["ibd_to_binary.py", "ibd_input.py"])
@pytest.mark.parametrize("pheno_mode", ["omitted", "excluded_and_missing", "nonexistent"])
@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("min_cm", [None, 3])
def test_conversion_keeps_all_samples(script, pheno_mode, compressed, min_cm, tmp_path):
    source = tmp_path / "segments.txt"
    content = (
        "case 1 control 1 1 100 500 1\n"
        "case 1 missing 1 1 200 400 4\n"
        "excluded 1 control 1 1 300 600 5\n"
        "missing 2 case 2 1 400 200 4\n"
    )
    source.write_bytes(gzip.compress(content.encode()) if compressed else content.encode())
    pheno = tmp_path / "pheno.txt"
    if pheno_mode == "excluded_and_missing":
        pheno.write_text("ID PHENO\ncase 1\ncontrol 0\nexcluded NA\n")
    # Use a path containing spaces to exercise text output without shell gzip.
    prefix = tmp_path / "converted data"
    command = [sys.executable, str(PROJECT_ROOT / "tools" / script),
               "-i", str(source), "-o", str(prefix), "-f", "hap-ibd", "-c", "1"]
    if pheno_mode != "omitted":
        command += ["--pheno", str(pheno)]
    if min_cm is not None:
        command += ["--min", str(min_cm)]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    if pheno_mode != "omitted":
        assert "--pheno is ignored" in result.stderr

    expected = {
        100: Counter({("case", "control"): 1}),
        200: Counter({("case", "control"): 1, ("case", "missing"): 2}),
        300: Counter({("case", "control"): 1, ("case", "missing"): 2, ("control", "excluded"): 1}),
        400: Counter({("case", "control"): 1, ("control", "excluded"): 1}),
        500: Counter({("control", "excluded"): 1}),
        600: Counter(),
    }
    if min_cm is not None:
        expected = {pos: counts - Counter({("case", "control"): 1})
                    for pos, counts in expected.items() if pos not in (100, 500)}
    actual = {}
    active = Counter()
    if script == "ibd_to_binary.py":
        names = Path(f"{prefix}.chr1.samples").read_text().splitlines()
        assert names == ["case", "control", "missing", "excluded"]
        with Path(f"{prefix}.chr1.ibdf").open("rb") as handle, mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as data:
            header = read_header(data)
            assert header.n_samples == 4
            for pos, offset, _ in read_index(data, header):
                kind, block = decode_block(data, offset)
                if kind == BLOCK_CHECKPOINT:
                    active = Counter(tuple(sorted((names[a], names[b]))) for a, b in zip(block[1], block[2]))
                else:
                    for a, b in zip(block[4], block[5]):
                        pair = tuple(sorted((names[a], names[b])))
                        assert active[pair] > 0
                        active[pair] -= 1
                        if not active[pair]:
                            del active[pair]
                    active.update(tuple(sorted((names[a], names[b]))) for a, b in zip(block[1], block[2]))
                actual[pos] = active.copy()
    else:
        with gzip.open(f"{prefix}.chr1.small.txt.gz", "rt") as handle:
            assert next(handle).strip() == "chr\tpos\tadd\tdel"
            for line in handle:
                _, pos, adds, dels = line.split("\t")
                for field, sign in [(adds, 1), (dels, -1)]:
                    for entry in field.split():
                        if entry == "NA":
                            continue
                        pair = tuple(sorted(entry.split(":", 1)[1].split("-")))
                        active[pair] += sign
                        if not active[pair]:
                            del active[pair]
                actual[int(pos)] = active.copy()
    assert actual == expected


@pytest.mark.parametrize("fmt,row", [
    ("hap-ibd", "a 1 b 1 1 100 200 1\n"),
    ("ilash", "a 1 b 1 1 100 200 x x 1\n"),
    ("germline", "a 1 b 1 1 100 200 x x x 1 cM\n"),
    ("rapid", "1 a b x x 100 200 1\n"),
    ("other:1;2;3;4;5;6", "a b 1 100 200 1\n"),
])
def test_sample_discovery_and_explicit_length_filter(fmt, row, tmp_path):
    source = tmp_path / "input.txt"
    source.write_text(row)
    segments, events, names = read_segments(str(source), Indices(fmt))
    assert names == ["a", "b"]
    assert segments == [(1.0, 0, 1)]
    assert set(events) == {100, 200}
    segments, events, names = read_segments(str(source), Indices(fmt), min_cm=3)
    assert segments == []
    assert not events
    assert names == ["a", "b"]


def test_binary_mappings_do_not_overwrite_other_chromosomes(tmp_path):
    prefix = tmp_path / "converted"
    for chrom, pair in [(1, ("a", "b")), (2, ("b", "c"))]:
        source = tmp_path / f"chr{chrom}.txt"
        source.write_text(f"{pair[0]} 1 {pair[1]} 1 {chrom} 100 200 4\n")
        subprocess.run([sys.executable, str(PROJECT_ROOT / "tools" / "ibd_to_binary.py"),
                        "-i", str(source), "-o", str(prefix), "-f", "hap-ibd", "-c", str(chrom)],
                       check=True, capture_output=True, text=True)
    assert Path(f"{prefix}.chr1.samples").read_text().splitlines() == ["a", "b"]
    assert Path(f"{prefix}.chr2.samples").read_text().splitlines() == ["b", "c"]
    assert not Path(f"{prefix}.samples").exists()
