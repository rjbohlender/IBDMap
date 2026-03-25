# IBDF v3 — IBD Binary Format Specification

**Version:** 3
**Date:** 2026-03-18
**Status:** Current
**Byte order:** Little-endian throughout
**Alignment:** 32 bytes (AVX2-compatible)

---

## 1. Overview

IBDF v3 is a binary format for storing Identity-by-Descent (IBD) pair data
at genomic breakpoints. It uses **checkpointed deltas**: most breakpoints
store only the segments added and removed since the previous breakpoint,
while every N-th breakpoint stores the complete set of active IBD pairs.

This design enables parallel processing. Any thread can load the nearest
checkpoint, replay a bounded number of deltas to reach its target position,
and then process its assigned range independently.

### 1.1 Design Goals

1. **Random access** — A position index at the end of the file allows
   binary search to any breakpoint.
2. **Bounded replay** — Checkpoint blocks cap the number of deltas any
   reader must replay before beginning useful work.
3. **SIMD-friendly** — Struct-of-arrays (SoA) layout within each block
   places all `cm`, `p1`, and `p2` values in contiguous, 32-byte-aligned
   arrays suitable for AVX2 / SSE4.1 vector loads.
4. **Zero-copy reads** — The format is designed for `mmap`; decoded views
   are pointers directly into mapped memory with no allocation or copying.

### 1.2 Companion Files

An IBDF file is accompanied by a plain-text `.samples` file that maps
the `uint32` sample indices used in the binary back to string sample
identifiers, one name per line, zero-indexed.

---

## 2. File Layout

```
Offset              Content
──────────────────────────────────────────────
0x00                File Header       (64 bytes, fixed)
0x40                Data Block 0      (variable, 32-byte aligned start)
...                 Data Block 1..N-1
(index_offset)      Position Index    (16 bytes × N)
EOF
```

The position index is written **after** all data blocks so the writer can
stream blocks without knowing final offsets in advance. Readers load the
index first (its offset is in the header), then seek to individual blocks.

---

## 3. File Header

**Size:** 64 bytes
**Offset:** 0

```
Offset  Size  Type      Field                 Description
──────────────────────────────────────────────────────────────────
0x00    4     uint32    magic                 0x33444249 ("IBD3" LE)
0x04    2     uint16    version               3
0x06    2     uint16    flags                 Reserved, must be 0
0x08    8     uint64    n_positions           Total number of breakpoints
0x10    8     uint64    index_offset          Byte offset to position index
0x18    4     uint32    n_samples             Number of unique sample IDs
0x1C    4     uint32    checkpoint_interval   Breakpoints between checkpoints
0x20    32    uint8[]   _reserved             Zero-filled, reserved for future use
```

**Total: 64 bytes**

**Validation:** Readers must verify `magic == 0x33444249` and
`version == 3` before proceeding. Unknown flags must be rejected.

---

## 4. Data Blocks

Each breakpoint is stored as one data block. Every block begins at a
**32-byte-aligned** file offset. If the previous block's data does not
end on a 32-byte boundary, zero padding is inserted.

### 4.1 Block Header

**Size:** 12 bytes

```
Offset  Size  Type      Field       Description
────────────────────────────────────────────────────────
+0x00   4     uint32    block_type  0 = DELTA, 1 = CHECKPOINT
+0x04   4     uint32    count_a     CHECKPOINT: n_pairs
                                    DELTA:      n_adds
+0x08   4     uint32    count_b     CHECKPOINT: 0 (unused)
                                    DELTA:      n_dels
```

### 4.2 Checkpoint Block (block_type = 1)

Contains the **complete set** of IBD pairs active at this breakpoint.
The first block in the file is always a checkpoint.

```
Offset              Size                    Content
──────────────────────────────────────────────────────────
+0x00               12                      BlockHeader { 1, n_pairs, 0 }
ALIGN_UP(+12)       n_pairs × 4            cm[n_pairs]     float32
ALIGN_UP(above)     n_pairs × 4            p1[n_pairs]     uint32
ALIGN_UP(above)     n_pairs × 4            p2[n_pairs]     uint32
```

Each array starts at a 32-byte-aligned offset. If `n_pairs` is 0
(no active IBD pairs), the arrays are empty but alignment padding
is still present between them.

### 4.3 Delta Block (block_type = 0)

Contains only the segments **added** and **removed** at this breakpoint
relative to the previous breakpoint's active set.

```
Offset              Size                    Content
──────────────────────────────────────────────────────────
+0x00               12                      BlockHeader { 0, n_adds, n_dels }
ALIGN_UP(+12)       n_adds × 4             add_cm[n_adds]     float32
ALIGN_UP(above)     n_adds × 4             add_p1[n_adds]     uint32
ALIGN_UP(above)     n_adds × 4             add_p2[n_adds]     uint32
ALIGN_UP(above)     n_dels × 4             del_cm[n_dels]     float32
ALIGN_UP(above)     n_dels × 4             del_p1[n_dels]     uint32
ALIGN_UP(above)     n_dels × 4             del_p2[n_dels]     uint32
```

**Add arrays** describe segments whose IBD sharing begins at this position.
**Del arrays** describe segments whose IBD sharing ends at this position.

### 4.4 Array Fields

Each IBD pair is represented by three parallel arrays in SoA layout:

| Array | Type      | Description                                        |
|-------|-----------|----------------------------------------------------|
| `cm`  | `float32` | centiMorgan length of the IBD segment              |
| `p1`  | `uint32`  | Numeric ID of the first sample (`p1 < p2` always)  |
| `p2`  | `uint32`  | Numeric ID of the second sample                    |

**Invariant:** `p1[i] < p2[i]` for all `i`. This canonical ordering
ensures each pair has a unique representation for matching during
delta removal.

**Precision:** cM values are stored as `float32` (IEEE 754 single
precision, ~7 significant digits). The original double-precision value
from the input file is narrowed during conversion.

### 4.5 Alignment Rules

The function `ALIGN_UP(x) = (x + 31) & ~31` rounds up to the next
32-byte boundary. It is applied:

1. **Before each block** — the block header starts at a 32-byte offset.
2. **Before each array** — every `cm`, `p1`, `p2` array within a block
   starts at a 32-byte offset from the beginning of the file.

If an array has 0 elements, no data is written but the offset still
advances to the next alignment boundary. This ensures subsequent arrays
retain their alignment guarantees.

**Padding bytes** are always zero-filled.

---

## 5. Position Index

**Location:** starts at byte offset `index_offset` (from the file header)
**Entry count:** `n_positions` (from the file header)
**Entry size:** 16 bytes each

```
Offset  Size  Type      Field        Description
─────────────────────────────────────────────────────────────
+0x00   8     uint64    bp_pos    Chromosome base-pair position
+0x08   8     uint64    data_offset  Byte offset to block (see below)
```

### 5.1 Checkpoint Flag

Bit 63 (the most significant bit) of `data_offset` is the **checkpoint
flag**:

```
  Bit 63 = 1  →  This position is a CHECKPOINT block
  Bit 63 = 0  →  This position is a DELTA block
```

The actual byte offset is obtained by masking off bit 63:

```
  real_offset = data_offset & 0x7FFFFFFFFFFFFFFF
```

This allows readers to scan the index for checkpoints without loading
block data, enabling fast lookup of the nearest prior checkpoint for
any target position.

### 5.2 Ordering

Index entries are stored in **ascending `bp_pos` order**, matching the
genomic order of breakpoints. This enables binary search by position.

### 5.3 First Entry

The first entry in the index is always a checkpoint. This guarantees
that replay can always find a starting point by scanning backward.

---

## 6. Reconstructing the Active Set

To obtain the full set of IBD pairs at position index `T`:

```
1. Scan backward from T in the position index to find the nearest
   entry where is_checkpoint(data_offset) is true. Call this index C.

2. Load the checkpoint block at C. This is the initial active set.

3. For each index i from C+1 to T inclusive:
     a. If block at i is a CHECKPOINT: replace the active set entirely.
     b. If block at i is a DELTA:
        - For each entry in del arrays: remove from active set.
        - For each entry in add arrays: insert into active set.

4. The active set now reflects position T.
```

**Worst-case replay:** `checkpoint_interval - 1` delta blocks.

---

## 7. Parallel Processing

The checkpoint + index design enables a straightforward parallel
decomposition:

```
1. mmap the file.
2. Read the position index (n_positions × 16 bytes at index_offset).
3. Assign each of K threads a contiguous range of position indices:
     Thread k processes indices [k×chunk, min((k+1)×chunk, n_positions)).
4. Each thread independently:
     a. Calls reconstruct_active_set(start_of_my_range).
     b. Steps forward through its range applying deltas.
     c. At each position, processes the active set (e.g., SIMD reduction
        over the cm array, filtering by pair, etc.).
```

No synchronization is needed between threads. Each thread touches only
its assigned data blocks plus at most one prior checkpoint.

---

## 8. Sample ID Mapping

The `.samples` companion file is a plain-text file with one sample name
per line:

```
SAMPLE_0000
SAMPLE_0001
SAMPLE_0002
...
```

Line number (0-indexed) corresponds to the `uint32` numeric ID used in
`p1` and `p2` arrays. This file is written by the converter and is not
embedded in the binary to keep the format simple and the binary
`mmap`-friendly.

---

## 9. Worked Example

A file with 5 breakpoints, `checkpoint_interval = 3`:

```
File offset   Block type     Content
──────────────────────────────────────────────────────
0x000         File Header    magic=0x33444249, n_positions=5,
                             checkpoint_interval=3, ...
0x040         CHECKPOINT     Pos 1000: 150 active pairs
0x4E0         DELTA          Pos 1050: +3 adds, -2 dels
0x560         DELTA          Pos 1100: +1 add,  -4 dels
0x5C0         CHECKPOINT     Pos 1200: 148 active pairs
0xA40         DELTA          Pos 1250: +2 adds, -1 del
0xAC0         Position Index 5 entries (80 bytes)
```

Position index contents:

| # | bp_pos | data_offset (hex)       | Checkpoint? |
|---|-----------|-------------------------|-------------|
| 0 | 1000      | `0x8000000000000040`    | Yes (bit 63 set) |
| 1 | 1050      | `0x00000000000004E0`    | No  |
| 2 | 1100      | `0x0000000000000560`    | No  |
| 3 | 1200      | `0x80000000000005C0`    | Yes (bit 63 set) |
| 4 | 1250      | `0x0000000000000A40`    | No  |

To reconstruct position 1100 (index 2): load checkpoint at index 0,
apply delta at index 1, apply delta at index 2.

To reconstruct position 1250 (index 4): load checkpoint at index 3,
apply delta at index 4. Only 1 replay step.

---

## 10. Size Characteristics

**Per-entry cost:** 12 bytes (4B `float32` + 4B `uint32` + 4B `uint32`).

**Block overhead per breakpoint:**

| Component              | Bytes                        |
|------------------------|------------------------------|
| Block header           | 12                           |
| Alignment to 32B       | 20 (header padding)          |
| Inter-array padding    | 0–31 per array (×3 or ×6)   |

**Checkpoint space cost** scales with `active_set_size × 12`.
**Delta space cost** scales with `(n_adds + n_dels) × 12`.

Typical overhead from checkpoints at interval N=1000 is 3–5% above
a pure-delta encoding.

---

## 11. Limitations and Constraints

- **Single chromosome per file.** Multi-chromosome data requires one
  IBDF file per chromosome.
- **Maximum 2^32 - 1 samples** (uint32 pair IDs).
- **Maximum 2^63 - 1 byte** file size (offset field is 63 bits after
  masking the checkpoint flag).
- **cM precision** is float32 (~7 significant digits). If higher
  precision is required, a future version could use float64 blocks.
- **No built-in compression.** Files can be compressed externally
  (e.g., with zstd) but lose `mmap` random-access capability.
  A future version could add per-block compression.
- **Pair identity for delta removal** relies on exact `(cm, p1, p2)`
  matching. If two segments share the same sample pair and identical
  float32 cM value, they are indistinguishable. In practice this is
  vanishingly rare but callers should be aware.

---

## 12. Reference Implementations

| Component    | File                | Language |
|--------------|---------------------|----------|
| Format spec  | `ibd_format.hpp`    | C++20    |
| Converter    | `ibd_to_binary.py`  | Python 3 |
| Reader       | `ibd_read.cpp`      | C++20    |
| Test suite   | `test_converter.py` | Python 3 |
