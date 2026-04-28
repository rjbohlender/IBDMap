#!/usr/bin/env python3
"""Convert legacy IBDMap text/zstd result output to Parquet."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import TextIO

import pyarrow as pa
import pyarrow.parquet as pq
import zstandard as zstd


SCHEMA = pa.schema(
    [
        pa.field("chrom", pa.string(), nullable=False),
        pa.field("pos", pa.int32(), nullable=False),
        pa.field("orig_cscs_rate", pa.float64(), nullable=False),
        pa.field("orig_cscn_rate", pa.float64(), nullable=False),
        pa.field("orig_cncn_rate", pa.float64(), nullable=False),
        pa.field("original", pa.float64(), nullable=False),
        pa.field("permutations", pa.list_(pa.float32()), nullable=False),
    ]
)


def compression_type(path: Path) -> str:
    with path.open("rb") as handle:
        magic = handle.read(4)
    if magic[:2] == b"\x1f\x8b":
        return "gzip"
    if magic == b"\x28\xb5\x2f\xfd":
        return "zstd"
    return "none"


def open_text(path: Path) -> TextIO:
    compressed = compression_type(path)
    if compressed == "gzip":
        return gzip.open(path, "rt")
    if compressed == "zstd":
        return zstd.open(path, "rt")
    return path.open("r", encoding="utf-8")


def normalize_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def empty_columns() -> dict[str, list]:
    return {
        "chrom": [],
        "pos": [],
        "orig_cscs_rate": [],
        "orig_cscn_rate": [],
        "orig_cncn_rate": [],
        "original": [],
        "permutations": [],
    }


def append_line(columns: dict[str, list], line: str, line_no: int) -> None:
    fields = line.strip().split()
    if not fields:
        return
    if len(fields) < 6:
        raise ValueError(f"line {line_no}: expected at least 6 columns, found {len(fields)}")

    columns["chrom"].append(normalize_chrom(fields[0]))
    columns["pos"].append(int(fields[1]))
    columns["orig_cscs_rate"].append(float(fields[2]))
    columns["orig_cscn_rate"].append(float(fields[3]))
    columns["orig_cncn_rate"].append(float(fields[4]))
    columns["original"].append(float(fields[5]))
    columns["permutations"].append([float(value) for value in fields[6:]])


def table_from_columns(columns: dict[str, list]) -> pa.Table:
    return pa.table(
        {
            "chrom": pa.array(columns["chrom"], type=pa.string()),
            "pos": pa.array(columns["pos"], type=pa.int32()),
            "orig_cscs_rate": pa.array(columns["orig_cscs_rate"], type=pa.float64()),
            "orig_cscn_rate": pa.array(columns["orig_cscn_rate"], type=pa.float64()),
            "orig_cncn_rate": pa.array(columns["orig_cncn_rate"], type=pa.float64()),
            "original": pa.array(columns["original"], type=pa.float64()),
            "permutations": pa.array(columns["permutations"], type=pa.list_(pa.float32())),
        },
        schema=SCHEMA,
    )


def convert(input_path: Path, output_path: Path, batch_size: int, compression: str) -> int:
    rows = 0
    columns = empty_columns()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open_text(input_path) as input_file, pq.ParquetWriter(
        output_path, SCHEMA, compression=compression
    ) as writer:
        for line_no, line in enumerate(input_file, start=1):
            append_line(columns, line, line_no)
            if len(columns["pos"]) >= batch_size:
                writer.write_table(table_from_columns(columns))
                rows += len(columns["pos"])
                columns = empty_columns()

        if columns["pos"]:
            writer.write_table(table_from_columns(columns))
            rows += len(columns["pos"])

    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert legacy IBDMap text/zstd result output to Parquet."
    )
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input text, .gz, or .zst result file.")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output .parquet file.")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=8192,
        help="Rows to buffer per Parquet write batch.",
    )
    parser.add_argument(
        "--compression",
        default="zstd",
        choices=["zstd", "snappy", "gzip", "brotli", "lz4", "none"],
        help="Parquet column compression codec.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    compression = None if args.compression == "none" else args.compression
    rows = convert(args.input, args.output, args.batch_size, compression)
    print(f"Wrote {rows} rows to {args.output}")


if __name__ == "__main__":
    main()
