import gzip
import mmap
from pathlib import Path
import subprocess
import sys
from collections import Counter

import pytest
import zstandard

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "tools"))

from compare_breakpoint_loading import decode_block, read_header, read_index
from ibd_text_to_binary import convert
from ibd_to_binary import BLOCK_CHECKPOINT


@pytest.fixture
def input_files(tmp_path):
    pheno = tmp_path / "pheno.txt"
    pheno.write_text("ID PHENO\na-one 1\nb-two 0\nc NA\nd -9\n")
    source = tmp_path / "input.txt"
    source.write_text(
        "chr\tpos\tadd\tdel\n"
        "1\t100\t4:a-one-c\t4:a-one-c\n"
        "chr1\t200\t4:a-one-b-two 4:a-one-b-two 2:a-one-c\tNA\n"
        "1\t250\t4:a-one-b-two 5:a-one-c\t4:b-two-a-one 5:a-one-c\n"
        "1\t300\t6:a-one-b-two\t4:a-one-b-two\n"
        "1\t350\tNA\tNA\n"
        "1\t400\tNA\t4:a-one-b-two 6:b-two-a-one 2:a-one-c\n"
        "1\t500\t4:a-one-c\t4:a-one-c\n"
    )
    return source, pheno


def read_states(path):
    states = {}
    active = Counter()
    with path.open("rb") as handle, mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as data:
        header = read_header(data)
        index = read_index(data, header)
        for pos, offset, is_checkpoint in index:
            kind, block = decode_block(data, offset)
            assert is_checkpoint == (kind == BLOCK_CHECKPOINT)
            if is_checkpoint:
                active = Counter(zip(*block))
            else:
                for segment in zip(*block[3:]):
                    assert active[segment] > 0, "IBDMap would reject this deletion"
                    active[segment] -= 1
                    if active[segment] == 0:
                        del active[segment]
                active.update(zip(*block[:3]))
            states[pos] = active.copy()
    return header, index, states


@pytest.mark.parametrize("compression", ["plain", "gzip", "zstd"])
@pytest.mark.parametrize("interval", [1, 2, 3, 1000])
@pytest.mark.parametrize("min_cm", [0, 3])
def test_conversion_preserves_breakpoints_and_active_segments(input_files, tmp_path, compression, interval, min_cm):
    source, pheno = input_files
    contents = source.read_bytes()
    if compression == "gzip":
        source.write_bytes(gzip.compress(contents))
    elif compression == "zstd":
        source.write_bytes(zstandard.ZstdCompressor().compress(contents))
    binary, samples, count = convert(source, tmp_path / "converted", "chr1", min_cm, interval, pheno)
    assert binary.name == "converted.chr1.ibdf"
    assert samples.name == "converted.chr1.samples"
    assert samples.read_text().splitlines() == ["a-one", "b-two", "c", "d"]
    assert count == 7
    header, index, states = read_states(binary)
    assert header.n_samples == 4
    assert header.checkpoint_interval == interval
    assert [checkpoint for _, _, checkpoint in index] == [i % interval == 0 for i in range(7)]
    short = Counter({(2.0, 0, 2): 1}) if min_cm == 0 else Counter()
    overlap = Counter({(4.0, 0, 1): 2}) + short
    replacement = Counter({(4.0, 0, 1): 1, (6.0, 0, 1): 1}) + short
    assert states == {
        100: Counter(), 200: overlap, 250: overlap, 300: replacement,
        350: replacement, 400: Counter(), 500: Counter(),
    }


@pytest.mark.parametrize("body,error", [
    ("", "no breakpoints"),
    ("1\t100\tbad-entry\tNA\n", "invalid segment entry"),
    ("1\t100\t4:unknown-c\tNA\n", "unknown sample pair"),
    ("1\t100\tNA\t4:a-one-c\n", "deletion without enough"),
    ("1\t100\tNA\tNA\n1\t99\tNA\tNA\n", "strictly increasing"),
    ("1\t100\tNA\tNA\n1\t100\tNA\tNA\n", "strictly increasing"),
    ("2\t100\tNA\tNA\n", "does not match"),
    ("1\t100\tnan:a-one-c\tNA\n", "invalid segment length"),
    ("1 100 NA NA\n", "tab-separated"),
])
def test_invalid_input_does_not_replace_existing_outputs(input_files, tmp_path, body, error):
    source, pheno = input_files
    source.write_text("chr\tpos\tadd\tdel\n" + body)
    binary, samples = tmp_path / "converted.chr1.ibdf", tmp_path / "converted.chr1.samples"
    binary.write_bytes(b"existing binary")
    samples.write_text("existing sample mapping")
    with pytest.raises(ValueError, match=error):
        convert(source, tmp_path / "converted", "1", pheno_path=pheno)
    assert binary.read_bytes() == b"existing binary"
    assert samples.read_text() == "existing sample mapping"
    assert not list(tmp_path.glob(".ibdf-*"))


def test_ambiguous_hyphenated_ids_fail(input_files, tmp_path):
    source, pheno = input_files
    pheno.write_text("ID PHENO\na 1\nb-c 0\na-b 1\nc 0\n")
    source.write_text("chr\tpos\tadd\tdel\n1\t100\t4:a-b-c\tNA\n")
    with pytest.raises(ValueError, match="ambiguous sample pair"):
        convert(source, tmp_path / "converted", "1", pheno_path=pheno)


@pytest.mark.parametrize("interval", [0, -1])
def test_invalid_checkpoint_interval(input_files, tmp_path, interval):
    source, pheno = input_files
    with pytest.raises(ValueError, match="checkpoint interval"):
        convert(source, tmp_path / "converted", "1", checkpoint_interval=interval, pheno_path=pheno)


@pytest.mark.parametrize("with_pheno", [False, True])
def test_command_line(input_files, tmp_path, with_pheno):
    source, pheno = input_files
    if not with_pheno:
        source.write_text(source.read_text().replace("a-one", "a").replace("b-two", "b"))
    command = [
        sys.executable, str(PROJECT_ROOT / "tools" / "ibd_text_to_binary.py"),
        "-i", str(source), "-o", str(tmp_path / "cli"),
        "-c", "1", "-N", "2",
    ]
    if with_pheno:
        command.extend(["-p", str(pheno)])
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "Wrote 7 breakpoints" in result.stdout
    assert read_states(tmp_path / "cli.chr1.ibdf")[0].n_positions == 7


@pytest.mark.parametrize("compression", ["plain", "gzip", "zstd"])
def test_inferred_names_match_explicit_names(input_files, tmp_path, compression):
    source, pheno = input_files
    # The first row introduces a and c; b is discovered after the first checkpoint.
    contents = source.read_text().replace("a-one", "a").replace("b-two", "b").encode()
    if compression == "gzip":
        contents = gzip.compress(contents)
    elif compression == "zstd":
        contents = zstandard.ZstdCompressor().compress(contents)
    source.write_bytes(contents)
    pheno.write_text("ID PHENO\na 1\nc NA\nb 0\n")
    inferred, inferred_samples, _ = convert(source, tmp_path / "inferred", "1", checkpoint_interval=2)
    explicit, explicit_samples, _ = convert(source, tmp_path / "explicit", "1", checkpoint_interval=2, pheno_path=pheno)
    assert inferred_samples.read_text().splitlines() == ["a", "c", "b"]
    assert inferred_samples.read_bytes() == explicit_samples.read_bytes()
    assert inferred.read_bytes() == explicit.read_bytes()
    assert read_states(inferred)[0].n_samples == 3


@pytest.mark.parametrize("pair,error", [
    ("a-b-c", "ambiguous sample pair"),
    ("abc", "invalid sample pair"),
    ("a-", "invalid sample pair"),
    ("-b", "invalid sample pair"),
])
def test_inference_rejects_ambiguous_or_invalid_pairs(tmp_path, pair, error):
    source = tmp_path / "input.txt"
    source.write_text("chr\tpos\tadd\tdel\n"
                      "1\t100\t4:a-b 4:b-c\tNA\n"
                      f"1\t200\t4:{pair}\tNA\n")
    # Earlier occurrences of a, b, and c do not establish where a-b-c splits.
    with pytest.raises(ValueError, match=f"input.txt:3: {error}.*--pheno"):
        convert(source, tmp_path / "failed", "1")
    assert not (tmp_path / "failed.chr1.ibdf").exists()
    assert not (tmp_path / "failed.chr1.samples").exists()


def test_chromosomes_keep_their_own_sample_mappings(tmp_path):
    first = tmp_path / "chr1.txt"
    first.write_text("chr\tpos\tadd\tdel\n"
                     "1\t100\t4:b-c\tNA\n1\t200\tNA\t4:b-c\n")
    second = tmp_path / "chr2.txt"
    second.write_text("chr\tpos\tadd\tdel\n"
                      "2\t100\t4:a-b 4:b-c\tNA\n2\t200\tNA\t4:a-b 4:b-c\n")
    prefix = tmp_path / "converted"
    binary1, samples1, _ = convert(first, prefix, "1")
    original_binary, original_samples = binary1.read_bytes(), samples1.read_bytes()
    binary2, samples2, _ = convert(second, prefix, "2")

    assert samples1.name == "converted.chr1.samples"
    assert samples2.name == "converted.chr2.samples"
    assert samples1.read_text().splitlines() == ["b", "c"]
    assert samples2.read_text().splitlines() == ["a", "b", "c"]
    assert samples1.read_bytes() == original_samples
    assert binary1.read_bytes() == original_binary
    assert not (tmp_path / "converted.samples").exists()

    # Shared pair b-c has different numeric IDs but still resolves correctly
    # when each chromosome uses its own companion mapping.
    for binary, samples, expected in [
        (binary1, samples1, {("b", "c")}),
        (binary2, samples2, {("a", "b"), ("b", "c")}),
    ]:
        names = samples.read_text().splitlines()
        header, _, states = read_states(binary)
        assert header.n_samples == len(names)
        assert {(names[p1], names[p2]) for _, p1, p2 in states[100]} == expected
        assert not states[200]
