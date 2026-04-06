#!/usr/bin/env python3
"""
ibd_to_binary.py — Convert IBD segment files to IBDF v3 checkpointed binary.

Every breakpoint stores binary DELTAS (which segments start/end).
Every N-th breakpoint is a CHECKPOINT with the full active IBD set.
Any downstream thread can start from the nearest checkpoint and replay
a small number of deltas, breaking the sequential dependency.

Supported input: GERMLINE, iLASH, RaPID, hap-ibd, or custom columns.

Usage:
    python3 ibd_to_binary.py -i segments.gz -p pheno.txt -o output \\
        -f hap-ibd -m 2.0 -c 22 [-N 1000]
"""

import os
import sys
import gzip
import struct
import argparse
from collections import defaultdict

# ── Binary format constants ──────────────────────────────────────────

IBDF3_MAGIC   = 0x33444249   # "IBD3" little-endian
IBDF3_VERSION = 3
ALIGN         = 32

BLOCK_DELTA      = 0
BLOCK_CHECKPOINT = 1

CHECKPOINT_FLAG = 1 << 63

def align_up(x: int) -> int:
    return (x + ALIGN - 1) & ~(ALIGN - 1)


# ── Index parsing ────────────────────────────────────────────────────

class Indices:
    """Column indices for various IBD segment file formats."""
    __slots__ = ('id1_indx', 'id2_indx', 'chr_indx', 'str_indx',
                 'end_indx', 'cM_indx', 'unit_indx')

    def __init__(self, fmt: str):
        self.unit_indx = None
        self.id1_indx = 0
        self.id2_indx = 2
        self.chr_indx = 4
        self.str_indx = 5
        self.end_indx = 6

        f = fmt.lower()
        if f == 'germline':
            self.cM_indx   = 10
            self.unit_indx = 11
        elif f == 'ilash':
            self.cM_indx = 9
        elif f in ('hap-ibd', 'hapibd'):
            self.cM_indx = 7
        elif f == 'rapid':
            self.chr_indx = 0
            self.id1_indx = 1
            self.id2_indx = 2
            self.cM_indx  = 7
        elif f.startswith('other:') or f.startswith('others:'):
            indx = f.split(':')[1].split(';')
            if len(indx) < 6:
                sys.exit('other format requires 6 indices: id1;id2;chr;start;end;cM')
            self.id1_indx = int(indx[0]) - 1
            self.id2_indx = int(indx[1]) - 1
            self.chr_indx = int(indx[2]) - 1
            self.str_indx = int(indx[3]) - 1
            self.end_indx = int(indx[4]) - 1
            self.cM_indx  = int(indx[5]) - 1
        else:
            sys.exit(f'Unrecognized format: {fmt}\n'
                     'Supported: GERMLINE / iLASH / RaPID / hap-ibd / '
                     'other:id1;id2;chr;start;end;cM')


# ── Phenotype / sample reading ───────────────────────────────────────

def read_phenotypes(pheno_path: str) -> tuple[dict[str, int], list[str]]:
    uniq_ids: dict[str, int] = {}
    id_list: list[str] = []
    duplicates = 0

    with open(pheno_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2 or parts[1] == 'NA':
                continue
            name = parts[0]
            if name in uniq_ids:
                duplicates += 1
            else:
                uniq_ids[name] = len(id_list)
                id_list.append(name)

    print(f'  Samples:       {len(id_list):,}')
    if duplicates:
        print(f'  Duplicates:    {duplicates:,} (skipped)')
    return uniq_ids, id_list


# ── Segment reading ──────────────────────────────────────────────────

def _is_gzipped(path: str) -> bool:
    with open(path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'


# A segment is stored as (cm, p1_idx, p2_idx) identified by its index
# in the segments list. Events map positions to lists of segment indices.

def read_segments(input_path: str, idx: Indices, uniq_ids: dict[str, int],
                  min_cm: float):
    """
    Returns:
        segments:  list of (cM, p1_idx, p2_idx)
        events:    {position: (add_seg_ids, remove_seg_ids)}
    """
    segments: list[tuple[float, int, int]] = []
    events: dict[int, tuple[list[int], list[int]]] = defaultdict(lambda: ([], []))

    skip_malformed  = 0
    skip_sample     = 0
    skip_min_cm     = 0
    skip_unit       = 0

    opener = gzip.open if _is_gzipped(input_path) else open
    with opener(input_path, 'rt') as f:
        for line in f:
            fields = line.split()
            try:
                id1   = fields[idx.id1_indx]
                id2   = fields[idx.id2_indx]
                cm    = float(fields[idx.cM_indx])
                start = int(fields[idx.str_indx])
                end   = int(fields[idx.end_indx])
            except (IndexError, ValueError):
                skip_malformed += 1
                continue

            if id1 not in uniq_ids or id2 not in uniq_ids:
                skip_sample += 1
                continue
            if cm < min_cm:
                skip_min_cm += 1
                continue
            if idx.unit_indx is not None and fields[idx.unit_indx] != 'cM':
                skip_unit += 1
                continue

            if start > end:
                start, end = end, start

            p1_idx = uniq_ids[id1]
            p2_idx = uniq_ids[id2]
            if p1_idx > p2_idx:
                p1_idx, p2_idx = p2_idx, p1_idx

            seg_id = len(segments)
            segments.append((cm, p1_idx, p2_idx))
            events[start][0].append(seg_id)
            events[end][1].append(seg_id)

    total_skipped = skip_malformed + skip_sample + skip_min_cm + skip_unit
    print(f'  Segments:      {len(segments):,}')
    if total_skipped:
        print(f'  Skipped:       {total_skipped:,}')
        if skip_malformed:
            print(f'    malformed:   {skip_malformed:,}')
        if skip_sample:
            print(f'    sample:      {skip_sample:,} (ID not in phenotype file)')
        if skip_min_cm:
            print(f'    min cM:      {skip_min_cm:,} (< {min_cm})')
        if skip_unit:
            print(f'    unit:        {skip_unit:,} (not cM)')
    print(f'  Breakpoints:   {len(events):,}')

    return segments, events


# ── Binary writer ────────────────────────────────────────────────────

class BinaryWriter:
    """Writes IBDF v3 checkpointed-delta binary."""

    # Header: magic(4) ver(2) flags(2) n_pos(8) idx_off(8) n_samp(4) ckpt_int(4) pad(32) = 64
    HEADER_FMT  = '<IHHQQII32s'
    HEADER_SIZE = 64

    INDEX_ENTRY_FMT  = '<QQ'
    INDEX_ENTRY_SIZE = 16

    BLOCK_HEADER_FMT  = '<III'   # type, count_a, count_b
    BLOCK_HEADER_SIZE = 12

    def __init__(self, path: str):
        self.path = path
        self.fout = open(path, 'wb')
        self.file_pos = 0
        self.index: list[tuple[int, int, bool]] = []  # (chrom_pos, offset, is_ckpt)

        # placeholder header
        hdr = struct.pack(self.HEADER_FMT,
                          IBDF3_MAGIC, IBDF3_VERSION, 0, 0, 0, 0, 0, b'\x00' * 32)
        self.fout.write(hdr)
        self.file_pos = self.HEADER_SIZE

    def _pad(self):
        aligned = align_up(self.file_pos)
        pad = aligned - self.file_pos
        if pad > 0:
            self.fout.write(b'\x00' * pad)
            self.file_pos = aligned

    def _write_array(self, values, fmt_char: str):
        """Write an aligned array of f32 or u32 values."""
        self._pad()
        n = len(values)
        if n > 0:
            buf = struct.pack(f'<{n}{fmt_char}', *values)
            self.fout.write(buf)
            self.file_pos += len(buf)

    def write_checkpoint(self, chrom_pos: int,
                         cm: list[float], p1: list[int], p2: list[int]):
        n = len(cm)
        self._pad()
        self.index.append((chrom_pos, self.file_pos, True))

        hdr = struct.pack(self.BLOCK_HEADER_FMT, BLOCK_CHECKPOINT, n, 0)
        self.fout.write(hdr)
        self.file_pos += self.BLOCK_HEADER_SIZE

        self._write_array(cm, 'f')
        self._write_array(p1, 'I')
        self._write_array(p2, 'I')

    def write_delta(self, chrom_pos: int,
                    add_cm: list[float], add_p1: list[int], add_p2: list[int],
                    del_cm: list[float], del_p1: list[int], del_p2: list[int]):
        na, nd = len(add_cm), len(del_cm)
        self._pad()
        self.index.append((chrom_pos, self.file_pos, False))

        hdr = struct.pack(self.BLOCK_HEADER_FMT, BLOCK_DELTA, na, nd)
        self.fout.write(hdr)
        self.file_pos += self.BLOCK_HEADER_SIZE

        self._write_array(add_cm, 'f')
        self._write_array(add_p1, 'I')
        self._write_array(add_p2, 'I')
        self._write_array(del_cm, 'f')
        self._write_array(del_p1, 'I')
        self._write_array(del_p2, 'I')

    def finalize(self, n_samples: int, checkpoint_interval: int):
        self._pad()
        index_offset = self.file_pos

        for chrom_pos, offset, is_ckpt in self.index:
            flagged_offset = offset | CHECKPOINT_FLAG if is_ckpt else offset
            buf = struct.pack(self.INDEX_ENTRY_FMT, chrom_pos, flagged_offset)
            self.fout.write(buf)

        n_positions = len(self.index)
        self.fout.seek(0)
        hdr = struct.pack(self.HEADER_FMT,
                          IBDF3_MAGIC, IBDF3_VERSION, 0,
                          n_positions, index_offset,
                          n_samples, checkpoint_interval, b'\x00' * 32)
        self.fout.write(hdr)
        self.fout.close()

        n_ckpt = sum(1 for _, _, c in self.index if c)
        return n_positions, index_offset, n_ckpt


# ── Sweep with checkpointed deltas ───────────────────────────────────

def sweep_and_write(segments, events, writer: BinaryWriter,
                    checkpoint_interval: int):
    """
    Sweep breakpoints in order. Every checkpoint_interval positions,
    write the full active set (checkpoint). Otherwise write deltas.
    """
    active: set[int] = set()
    sorted_positions = sorted(events.keys())

    max_active = 0
    total_delta_entries = 0
    total_checkpoint_entries = 0
    positions_since_checkpoint = checkpoint_interval  # force first to be checkpoint

    for pos in sorted_positions:
        adds_ids, removes_ids = events[pos]

        # apply removals
        for seg_id in removes_ids:
            active.discard(seg_id)
        # apply additions
        for seg_id in adds_ids:
            active.add(seg_id)

        n_active = len(active)
        if n_active > max_active:
            max_active = n_active

        if positions_since_checkpoint >= checkpoint_interval:
            # ── Write checkpoint: full active set ────────────────────
            cm_arr  = [segments[s][0] for s in active]
            p1_arr  = [segments[s][1] for s in active]
            p2_arr  = [segments[s][2] for s in active]
            writer.write_checkpoint(pos, cm_arr, p1_arr, p2_arr)
            total_checkpoint_entries += n_active
            positions_since_checkpoint = 1
        else:
            # ── Write delta: just adds and dels ──────────────────────
            add_cm  = [segments[s][0] for s in adds_ids]
            add_p1  = [segments[s][1] for s in adds_ids]
            add_p2  = [segments[s][2] for s in adds_ids]
            del_cm  = [segments[s][0] for s in removes_ids]
            del_p1  = [segments[s][1] for s in removes_ids]
            del_p2  = [segments[s][2] for s in removes_ids]
            writer.write_delta(pos, add_cm, add_p1, add_p2,
                               del_cm, del_p1, del_p2)
            total_delta_entries += len(adds_ids) + len(removes_ids)
            positions_since_checkpoint += 1

    return max_active, total_delta_entries, total_checkpoint_entries


# ── Main ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Convert IBD segment files to IBDF v3 checkpointed binary.')
    parser.add_argument('-i', '--input', required=True,
                        help='Input IBD segment file (plain or gzipped).')
    parser.add_argument('-p', '--pheno', required=True,
                        help='Phenotype file (sample filtering + ID assignment).')
    parser.add_argument('-o', '--output', required=True,
                        help='Output prefix → <prefix>.chr<c>.ibdf + <prefix>.samples')
    parser.add_argument('-f', '--format', required=True,
                        help='Input format: GERMLINE / iLASH / RaPID / hap-ibd / '
                             'other:id1;id2;chr;start;end;cM')
    parser.add_argument('-m', '--min', type=float, required=True,
                        help='Minimum cM length to include.')
    parser.add_argument('-c', '--chrom', required=True,
                        help='Chromosome label.')
    parser.add_argument('-N', '--checkpoint-interval', type=int, default=1000,
                        help='Write a full checkpoint every N breakpoints '
                             '(default: 1000). Lower = faster thread startup, '
                             'larger file. Higher = smaller file, more replay.')
    args = parser.parse_args()

    ckpt_n = args.checkpoint_interval

    print(f'IBD to Binary Converter (IBDF v{IBDF3_VERSION}, checkpointed)')
    print(f'  Input:         {args.input}')
    print(f'  Phenotype:     {args.pheno}')
    print(f'  Output:        {args.output}.chr{args.chrom}.ibdf')
    print(f'  Format:        {args.format}')
    print(f'  Min cM:        {args.min}')
    print(f'  Chromosome:    {args.chrom}')
    print(f'  Checkpoint:    every {ckpt_n} breakpoints')
    print()

    idx = Indices(args.format)

    print('Reading phenotypes...')
    uniq_ids, id_list = read_phenotypes(args.pheno)

    print('Reading IBD segments...')
    segments, events = read_segments(args.input, idx, uniq_ids, args.min)

    if not segments:
        print('No segments passed filters. Exiting.')
        sys.exit(0)

    bin_path = f'{args.output}.chr{args.chrom}.ibdf'
    print(f'Writing binary to {bin_path}...')

    writer = BinaryWriter(bin_path)
    max_active, total_delta, total_ckpt = sweep_and_write(
        segments, events, writer, ckpt_n)
    n_positions, index_offset, n_ckpt = writer.finalize(
        n_samples=len(id_list), checkpoint_interval=ckpt_n)

    # sample mapping file
    samples_path = f'{args.output}.samples'
    with open(samples_path, 'w') as f:
        for name in id_list:
            f.write(f'{name}\n')

    bin_size = os.path.getsize(bin_path)
    n_deltas = n_positions - n_ckpt
    total_entries = total_delta + total_ckpt

    print()
    print(f'  Breakpoints:   {n_positions:,}')
    print(f'  Checkpoints:   {n_ckpt:,} (every {ckpt_n})')
    print(f'  Deltas:        {n_deltas:,}')
    print(f'  Max active:    {max_active:,} pairs')
    print(f'  Delta entries: {total_delta:,}')
    print(f'  Ckpt entries:  {total_ckpt:,}')
    print(f'  Total entries: {total_entries:,}')
    print(f'  Binary size:   {bin_size / 1e6:.1f} MB')
    print(f'  Samples file:  {samples_path}')
    print()

    # compare: what would pure delta (v1) and pure full-set (v2) cost?
    # v1: every position is a delta → total_delta + total_ckpt is wrong;
    # just use 2 * n_segments (each segment has exactly 1 add + 1 del)
    pure_delta_entries = 2 * len(segments)
    pure_delta_bytes = pure_delta_entries * 12 + n_positions * (96 + 16)  # rough
    print(f'  vs pure delta: ~{pure_delta_entries:,} entries '
          f'(~{pure_delta_bytes / 1e6:.1f} MB)')
    print(f'  vs full-set:   would be ~{total_ckpt / max(n_ckpt,1) * n_positions:,.0f} '
          f'entries')
    print()
    print('Done.')


if __name__ == '__main__':
    main()
