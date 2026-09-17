#!/usr/bin/env python3
"""Convert legacy IBDMap chr/pos/add/del text input to IBDF v3."""

import argparse
from collections import Counter
import gzip
import math
import os
from pathlib import Path
import tempfile
from typing import Optional

from ibd_to_binary import BinaryWriter


def open_text(path: Path):
    with path.open("rb") as handle:
        magic = handle.read(4)
    if magic[:2] == b"\x1f\x8b":
        return gzip.open(path, "rt")
    if magic == b"\x28\xb5\x2f\xfd":
        try:
            import zstandard
        except ImportError as exc:
            raise RuntimeError("Zstandard input requires `pip install zstandard`.") from exc
        return zstandard.open(path, "rt")
    return path.open(encoding="utf-8")


def read_sample_names(path: Path) -> dict[str, int]:
    """Read IDs from a phenotype file with a header; retain all phenotype statuses."""
    names = {}
    header_seen = False
    with open_text(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if not header_seen:
                header_seen = True
                continue
            name = line.split()[0]
            if name in names:
                raise ValueError(f"{path}:{lineno}: duplicate sample ID {name!r}")
            names[name] = len(names)
    if not names:
        raise ValueError(f"{path}: no sample IDs found after the phenotype header")
    return names


def parse_events(field: str, names: dict[str, int], min_cm: float,
                 infer_names: bool = False) -> Counter:
    result = Counter()
    for entry in field.split():
        if entry == "NA":
            if field.strip() != "NA":
                raise ValueError("NA must appear alone in an event column")
            continue
        try:
            length, pair = entry.split(":", 1)
            cm = float(length)
        except ValueError as exc:
            raise ValueError(f"invalid segment entry {entry!r}; expected cM:sample1-sample2") from exc
        if not math.isfinite(cm) or cm < 0:
            raise ValueError(f"invalid segment length in {entry!r}")
        # Without a supplied sample list, every nonempty split is possible.
        # Do not use names inferred from earlier rows to guess an ambiguous pair.
        matches = [
            (pair[:idx], pair[idx + 1:])
            for idx, char in enumerate(pair)
            if char == "-" and 0 < idx < len(pair) - 1
        ]
        if not infer_names:
            matches = [(first, second) for first, second in matches
                       if first in names and second in names]
        if len(matches) != 1:
            reason = "ambiguous" if matches else ("invalid" if infer_names else "unknown")
            hint = "provide --pheno to resolve sample IDs" if infer_names else "check the phenotype sample IDs"
            raise ValueError(f"{reason} sample pair {pair!r}; {hint}")
        first, second = matches[0]
        if infer_names:
            for name in (first, second):
                if name not in names:
                    names[name] = len(names)
        if cm >= min_cm:
            p1, p2 = sorted((names[first], names[second]))
            result[(cm, p1, p2)] += 1
    return result


def segment_arrays(counts: Counter) -> tuple[list, list, list]:
    cm, p1, p2 = [], [], []
    for (length, first, second), count in counts.items():
        cm.extend([length] * count)
        p1.extend([first] * count)
        p2.extend([second] * count)
    return cm, p1, p2


def write_breakpoints(input_path: Path, writer: BinaryWriter, names: dict[str, int],
                      chrom: str, min_cm: float, checkpoint_interval: int,
                      infer_names: bool = False) -> int:
    active = Counter()
    previous_pos = -1
    header_seen = False
    n_positions = 0
    with open_text(input_path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if not header_seen:
                if [v.lower() for v in line.split()] != ["chr", "pos", "add", "del"]:
                    raise ValueError(f"{input_path}:{lineno}: expected chr pos add del header")
                header_seen = True
                continue
            try:
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 4:
                    raise ValueError("expected four tab-separated columns: chr, pos, add, del")
                row_chrom, position, add_field, del_field = fields
                if row_chrom.removeprefix("chr") != chrom.removeprefix("chr"):
                    raise ValueError(f"chromosome {row_chrom!r} does not match --chrom {chrom!r}")
                pos = int(position)
                if not 0 <= pos < 2**64 or pos <= previous_pos:
                    raise ValueError("positions must be nonnegative, unique, and strictly increasing")
                adds = parse_events(add_field, names, min_cm, infer_names)
                dels = parse_events(del_field, names, min_cm, infer_names)
                # Equal additions/deletions have no net effect. Cancel them
                # before writing deltas so readers never subtract a new segment.
                cancelled = adds & dels
                adds -= cancelled
                dels -= cancelled
                for segment, count in dels.items():
                    if active[segment] < count:
                        raise ValueError(f"deletion without enough active segments: {segment}")
                    active[segment] -= count
                    if not active[segment]:
                        del active[segment]
                active.update(adds)
                if n_positions % checkpoint_interval == 0:
                    writer.write_checkpoint(pos, *segment_arrays(active))
                else:
                    writer.write_delta(pos, *segment_arrays(adds), *segment_arrays(dels))
                previous_pos = pos
                n_positions += 1
            except (ValueError, OverflowError) as exc:
                raise ValueError(f"{input_path}:{lineno}: {exc}") from exc
    if not n_positions:
        raise ValueError(f"{input_path}: no breakpoints found")
    return n_positions


def convert(input_path: Path, output_prefix: Path, chrom: str,
            min_cm: float = 0.0, checkpoint_interval: int = 1000,
            pheno_path: Optional[Path] = None) -> tuple[Path, Path, int]:
    if checkpoint_interval < 1:
        raise ValueError("checkpoint interval must be at least 1")
    if not math.isfinite(min_cm) or min_cm < 0:
        raise ValueError("minimum cM must be finite and nonnegative")
    names = read_sample_names(pheno_path) if pheno_path is not None else {}
    chrom = chrom.removeprefix("chr")
    if not chrom or any(c in chrom for c in "/\\\t\n"):
        raise ValueError("invalid chromosome label")
    binary_path = Path(f"{output_prefix}.chr{chrom}.ibdf")
    # Inferred IDs depend on which samples appear on this chromosome. Keep
    # each mapping beside its binary file so another conversion cannot replace it.
    samples_path = Path(f"{output_prefix}.chr{chrom}.samples")
    binary_path.parent.mkdir(parents=True, exist_ok=True)
    # Publish only after successful validation, leaving existing outputs intact
    # if a malformed row occurs late in a large input file.
    with tempfile.TemporaryDirectory(prefix=".ibdf-", dir=binary_path.parent) as tmp:
        staged_binary = Path(tmp) / "output.ibdf"
        staged_samples = Path(tmp) / "output.samples"
        writer = BinaryWriter(str(staged_binary))
        try:
            n_positions = write_breakpoints(input_path, writer, names, chrom, min_cm,
                                            checkpoint_interval, infer_names=pheno_path is None)
            writer.finalize(len(names), checkpoint_interval)
        finally:
            writer.fout.close()
        staged_samples.write_text("".join(f"{name}\n" for name in names), encoding="utf-8")
        os.replace(staged_binary, binary_path)
        os.replace(staged_samples, samples_path)
    return binary_path, samples_path, n_positions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, type=Path,
                        help="IBDMap chr/pos/add/del text file (plain, gzip, or Zstandard).")
    parser.add_argument("-p", "--pheno", type=Path,
                        help="Optional phenotype file to resolve ambiguous sample pairs (e.g. hyphenated IDs). "
                             "Requires a header and sample IDs in the first column; all statuses are retained. "
                             "Without it, sample IDs are inferred from unambiguous pairs.")
    parser.add_argument("-o", "--output", required=True, type=Path,
                        help="Output prefix for <prefix>.chr<c>.ibdf and <prefix>.chr<c>.samples.")
    parser.add_argument("-c", "--chrom", required=True, help="Chromosome label, with or without chr.")
    parser.add_argument("-m", "--min", type=float, default=0.0,
                        help="Minimum segment length in cM (default: retain all lengths).")
    parser.add_argument("-N", "--checkpoint-interval", type=int, default=1000)
    args = parser.parse_args()
    try:
        binary_path, samples_path, n_positions = convert(
            args.input, args.output, args.chrom, args.min, args.checkpoint_interval, args.pheno)
    except (ValueError, RuntimeError, OSError) as exc:
        parser.exit(1, f"Error: {exc}\n")
    print(f"Wrote {n_positions} breakpoints to {binary_path}\nSample mapping: {samples_path}")


if __name__ == "__main__":
    main()
