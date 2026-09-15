#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
import gzip
import mmap
import struct
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


MAGIC = 0x33444249
CHECKPOINT_FLAG = 1 << 63
BLOCK_DELTA = 0
BLOCK_CHECKPOINT = 1
ALIGN = 32


def align_up(value: int) -> int:
    return (value + ALIGN - 1) & ~(ALIGN - 1)


def read_phenotypes(path: Path, pheno_col: int = 1, swap: bool = False) -> tuple[list[str], list[int], dict[str, int]]:
    samples: list[str] = []
    phenotypes: list[int] = []
    lookup: dict[str, int] = {}
    seen: set[str] = set()

    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle):
            if line.startswith("#") or line_no == 0:
                continue
            fields = line.split()
            if len(fields) <= pheno_col:
                continue
            if fields[pheno_col] in {"NA", "-9"}:
                continue
            sample = fields[0]
            if sample in seen:
                raise ValueError(f"duplicate sample ID: {sample}")
            value = int(fields[pheno_col])
            if value not in {0, 1}:
                raise ValueError(f"invalid phenotype on line {line_no + 1}: {fields[pheno_col]}")
            if swap:
                value = 1 - value
            seen.add(sample)
            samples.append(sample)
            phenotypes.append(value)
            lookup[sample] = value

    order = sorted(range(len(samples)), key=samples.__getitem__)
    sorted_samples = [samples[i] for i in order]
    sorted_phenotypes = [phenotypes[i] for i in order]
    sorted_lookup = dict(zip(sorted_samples, sorted_phenotypes))
    return sorted_samples, sorted_phenotypes, sorted_lookup


def translate(samples: list[str], a: str, b: str) -> int:
    if a == b:
        return -1
    ia = bisect.bisect_left(samples, a)
    ib = bisect.bisect_left(samples, b)
    if ia == len(samples) or ib == len(samples) or samples[ia] != a or samples[ib] != b:
        return -1
    if ia > ib:
        ia, ib = ib, ia
    n = len(samples)
    start = (2 * n - 1 - ia) * ia // 2
    return start + ib - ia - 1


def back_translate(n: int, row: int) -> tuple[int, int]:
    # Inverse of IBDMap's row-major upper-triangle mapping.
    lo = 0
    hi = n - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        start = (2 * n - 1 - mid) * mid // 2
        if start <= row:
            lo = mid
        else:
            hi = mid - 1
    i = lo
    j = row - (2 * n - 1 - i) * i // 2 + i + 1
    return i, j


def counts_for_rows(rows: set[int], samples: list[str], phenotypes: list[int]) -> tuple[int, int, int]:
    cscs = cscn = cncn = 0
    n = len(samples)
    for row in rows:
        i, j = back_translate(n, row)
        x = phenotypes[i]
        y = phenotypes[j]
        cscs += x & y
        cscn += x ^ y
        cncn += (1 - x) & (1 - y)
    return cscs, cscn, cncn


def apply_text_entries(active: Counter[int], samples: list[str], field: str, sign: int, min_cm: float) -> tuple[int, int]:
    if field == "NA" or not field:
        return 0, 0
    applied = skipped = 0
    for entry in field.split():
        if entry == "NA":
            break
        try:
            cm_text, pair_text = entry.split(":", 1)
            cm = float(cm_text)
            a, b = pair_text.split("-", 1)
        except ValueError:
            skipped += 1
            continue
        if cm < min_cm:
            skipped += 1
            continue
        row = translate(samples, a, b)
        if row < 0:
            skipped += 1
            continue
        active[row] += sign
        if active[row] == 0:
            del active[row]
        applied += 1
    return applied, skipped


@dataclass
class TextResult:
    rows: set[int]
    counter: Counter[int]
    line_no: int
    pos: int
    add_applied: int
    del_applied: int
    add_skipped: int
    del_skipped: int


def load_text_target(path: Path, samples: list[str], target_pos: int, min_cm: float) -> TextResult:
    active: Counter[int] = Counter()
    with gzip.open(path, "rt") as handle:
        header = next(handle)
        if not header.startswith("chr"):
            raise ValueError(f"unexpected text header: {header[:80]}")
        for line_no, line in enumerate(handle, start=2):
            chrom, pos_text, add_field, del_field = line.rstrip("\n").split("\t", 3)
            pos = int(pos_text)
            add_applied, add_skipped = apply_text_entries(active, samples, add_field, 1, min_cm)
            del_applied, del_skipped = apply_text_entries(active, samples, del_field, -1, min_cm)
            if pos == target_pos:
                return TextResult(
                    rows={row for row, count in active.items() if count > 0},
                    counter=Counter({row: count for row, count in active.items() if count > 0}),
                    line_no=line_no,
                    pos=pos,
                    add_applied=add_applied,
                    del_applied=del_applied,
                    add_skipped=add_skipped,
                    del_skipped=del_skipped,
                )
            if pos > target_pos:
                raise ValueError(f"text stream passed target at {pos} without finding {target_pos}")
    raise ValueError(f"target {target_pos} not found in text stream")


@dataclass
class BinaryHeader:
    n_positions: int
    index_offset: int
    n_samples: int
    checkpoint_interval: int


def read_header(mm: mmap.mmap) -> BinaryHeader:
    magic, version, _flags, n_positions, index_offset, n_samples, checkpoint_interval, _pad = struct.unpack_from(
        "<IHHQQII32s", mm, 0
    )
    if magic != MAGIC or version != 3:
        raise ValueError(f"unexpected IBDF header magic={magic:x} version={version}")
    return BinaryHeader(n_positions, index_offset, n_samples, checkpoint_interval)


def read_index(mm: mmap.mmap, header: BinaryHeader) -> list[tuple[int, int, bool]]:
    out: list[tuple[int, int, bool]] = []
    for i in range(header.n_positions):
        pos, raw_offset = struct.unpack_from("<QQ", mm, header.index_offset + i * 16)
        is_checkpoint = (raw_offset & CHECKPOINT_FLAG) != 0
        out.append((pos, raw_offset & ~CHECKPOINT_FLAG, is_checkpoint))
    return out


def read_array(mm: mmap.mmap, offset: int, count: int, fmt: str) -> tuple[list[float] | list[int], int]:
    offset = align_up(offset)
    if count == 0:
        return [], offset
    size = 4 * count
    values = list(struct.unpack_from(f"<{count}{fmt}", mm, offset))
    return values, offset + align_up(size)


def decode_block(mm: mmap.mmap, offset: int):
    block_type, count_a, count_b = struct.unpack_from("<III", mm, offset)
    cursor = offset + align_up(12)
    if block_type == BLOCK_CHECKPOINT:
        cm, cursor = read_array(mm, cursor, count_a, "f")
        p1, cursor = read_array(mm, cursor, count_a, "I")
        p2, cursor = read_array(mm, cursor, count_a, "I")
        return block_type, (cm, p1, p2)
    if block_type == BLOCK_DELTA:
        add_cm, cursor = read_array(mm, cursor, count_a, "f")
        add_p1, cursor = read_array(mm, cursor, count_a, "I")
        add_p2, cursor = read_array(mm, cursor, count_a, "I")
        del_cm, cursor = read_array(mm, cursor, count_b, "f")
        del_p1, cursor = read_array(mm, cursor, count_b, "I")
        del_p2, cursor = read_array(mm, cursor, count_b, "I")
        return block_type, (add_cm, add_p1, add_p2, del_cm, del_p1, del_p2)
    raise ValueError(f"unknown block type {block_type}")


@dataclass
class BinaryResult:
    rows: set[int]
    counter: Counter[int]
    index_no: int
    pos: int
    is_checkpoint: bool
    adds: int
    dels: int
    add_skipped: int
    del_skipped: int
    checkpoint_size: int


def load_binary_target(
    ibdf_path: Path,
    binary_samples_path: Path,
    pheno_samples: list[str],
    target_pos: int,
    min_cm: float,
) -> BinaryResult:
    binary_samples = [line.strip() for line in binary_samples_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    with ibdf_path.open("rb") as handle:
        mm = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            header = read_header(mm)
            if header.n_samples != len(binary_samples):
                raise ValueError(f"IBDF says {header.n_samples} samples but sample file has {len(binary_samples)}")
            index = read_index(mm, header)
            positions = [pos for pos, _offset, _is_checkpoint in index]
            idx = bisect.bisect_left(positions, target_pos)
            if idx == len(index) or positions[idx] != target_pos:
                raise ValueError(f"target {target_pos} not found in IBDF index")

            chunk_start = idx
            while chunk_start > 0 and not index[chunk_start][2]:
                chunk_start -= 1
            if not index[chunk_start][2]:
                raise ValueError("no checkpoint before target")

            active: Counter[int] = Counter()
            last_adds = last_dels = last_add_skipped = last_del_skipped = 0
            checkpoint_size = 0
            for block_idx in range(chunk_start, idx + 1):
                pos, offset, is_checkpoint = index[block_idx]
                block_type, block = decode_block(mm, offset)
                if block_type == BLOCK_CHECKPOINT:
                    active.clear()
                    cm, p1, p2 = block
                    checkpoint_size = len(cm)
                    last_adds = last_dels = last_add_skipped = last_del_skipped = 0
                    for seg_cm, i1, i2 in zip(cm, p1, p2):
                        if seg_cm < min_cm:
                            last_add_skipped += 1
                            continue
                        row = translate(pheno_samples, binary_samples[i1], binary_samples[i2])
                        if row < 0:
                            last_add_skipped += 1
                            continue
                        active[row] += 1
                        last_adds += 1
                else:
                    add_cm, add_p1, add_p2, del_cm, del_p1, del_p2 = block
                    last_adds = last_dels = last_add_skipped = last_del_skipped = 0
                    # Current parser applies deletions before additions.
                    for seg_cm, i1, i2 in zip(del_cm, del_p1, del_p2):
                        if seg_cm < min_cm:
                            last_del_skipped += 1
                            continue
                        row = translate(pheno_samples, binary_samples[i1], binary_samples[i2])
                        if row < 0:
                            last_del_skipped += 1
                            continue
                        active[row] -= 1
                        if active[row] == 0:
                            del active[row]
                        last_dels += 1
                    for seg_cm, i1, i2 in zip(add_cm, add_p1, add_p2):
                        if seg_cm < min_cm:
                            last_add_skipped += 1
                            continue
                        row = translate(pheno_samples, binary_samples[i1], binary_samples[i2])
                        if row < 0:
                            last_add_skipped += 1
                            continue
                        active[row] += 1
                        last_adds += 1

            return BinaryResult(
                rows={row for row, count in active.items() if count > 0},
                counter=Counter({row: count for row, count in active.items() if count > 0}),
                index_no=idx,
                pos=target_pos,
                is_checkpoint=index[idx][2],
                adds=last_adds,
                dels=last_dels,
                add_skipped=last_add_skipped,
                del_skipped=last_del_skipped,
                checkpoint_size=checkpoint_size if index[idx][2] else 0,
            )
        finally:
            mm.close()


def describe(label: str, rows: set[int], counter: Counter[int], samples: list[str], phenotypes: list[int]) -> None:
    cscs, cscn, cncn = counts_for_rows(rows, samples, phenotypes)
    total_segment_count = sum(counter.values())
    duplicate_rows = sum(1 for value in counter.values() if value > 1)
    max_multiplicity = max(counter.values(), default=0)
    print(f"{label}:")
    print(f"  active nonzero sample-pair rows: {len(rows)}")
    print(f"  active segment multiplicity sum: {total_segment_count}")
    print(f"  rows with multiplicity >1: {duplicate_rows}")
    print(f"  max row multiplicity: {max_multiplicity}")
    print(f"  cscs/cscn/cncn over nonzero rows: {cscs}/{cscn}/{cncn}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibdf", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--text", required=True, type=Path)
    parser.add_argument("--pheno", required=True, type=Path)
    parser.add_argument("--target", required=True, type=int)
    parser.add_argument("--pheno-col", type=int, default=1)
    parser.add_argument("--min-cm", type=float, default=3.0)
    args = parser.parse_args()

    samples, phenotypes, _lookup = read_phenotypes(args.pheno, args.pheno_col)
    n_cases = sum(phenotypes)
    n_controls = len(phenotypes) - n_cases
    print(f"Phenotype samples after filtering/sorting: {len(samples)}")
    print(f"Cases/controls: {n_cases}/{n_controls}")

    text = load_text_target(args.text, samples, args.target, args.min_cm)
    binary = load_binary_target(args.ibdf, args.samples, samples, args.target, args.min_cm)

    print()
    print(f"Text target line: {text.line_no}, pos={text.pos}")
    print(
        f"  target-line applied add/del: {text.add_applied}/{text.del_applied}; "
        f"skipped add/del: {text.add_skipped}/{text.del_skipped}"
    )
    describe("Text active state", text.rows, text.counter, samples, phenotypes)

    print()
    print(
        f"Binary target index: {binary.index_no}, pos={binary.pos}, "
        f"is_checkpoint={binary.is_checkpoint}"
    )
    if binary.is_checkpoint:
        print(f"  checkpoint raw segment count: {binary.checkpoint_size}")
    print(
        f"  target-block applied add/del: {binary.adds}/{binary.dels}; "
        f"skipped add/del: {binary.add_skipped}/{binary.del_skipped}"
    )
    describe("Binary active state", binary.rows, binary.counter, samples, phenotypes)

    print()
    only_text = text.rows - binary.rows
    only_binary = binary.rows - text.rows
    print(f"Rows only in text: {len(only_text)}")
    print(f"Rows only in binary: {len(only_binary)}")
    if only_text or only_binary:
        print("First rows only in text:")
        for row in sorted(only_text)[:10]:
            i, j = back_translate(len(samples), row)
            print(f"  {row}: {samples[i]}-{samples[j]} phenos={phenotypes[i]}-{phenotypes[j]} text_count={text.counter[row]}")
        print("First rows only in binary:")
        for row in sorted(only_binary)[:10]:
            i, j = back_translate(len(samples), row)
            print(f"  {row}: {samples[i]}-{samples[j]} phenos={phenotypes[i]}-{phenotypes[j]} binary_count={binary.counter[row]}")

    count_mismatches = [
        row for row in sorted(text.rows & binary.rows)
        if text.counter[row] != binary.counter[row]
    ]
    print(f"Shared rows with different multiplicity: {len(count_mismatches)}")
    for row in count_mismatches[:10]:
        i, j = back_translate(len(samples), row)
        print(
            f"  {row}: {samples[i]}-{samples[j]} phenos={phenotypes[i]}-{phenotypes[j]} "
            f"text_count={text.counter[row]} binary_count={binary.counter[row]}"
        )


if __name__ == "__main__":
    main()
