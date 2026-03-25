#pragma once

#include <cstdint>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <vector>

/*
 * IBD Binary Format (IBDF) v3 — Checkpointed Deltas
 *
 * Every breakpoint stores binary DELTAS (adds/dels).
 * Every N-th breakpoint is a CHECKPOINT with the full active set.
 * Any thread can start from the nearest checkpoint and replay forward.
 *
 * File layout:
 *   [File Header]       64 bytes
 *   [Data Blocks]       one per breakpoint, 32-byte aligned
 *   [Position Index]    array of index entries at end of file
 *
 * Block types:
 *   CHECKPOINT (full active set):
 *     BlockHeader { type=1, n_pairs, 0 }
 *     cm[n_pairs]  f32, 32B aligned
 *     p1[n_pairs]  u32, 32B aligned
 *     p2[n_pairs]  u32, 32B aligned
 *
 *   DELTA (adds + dels since previous breakpoint):
 *     BlockHeader { type=0, n_adds, n_dels }
 *     add_cm[n_adds]  f32    del_cm[n_dels]  f32
 *     add_p1[n_adds]  u32    del_p1[n_dels]  u32
 *     add_p2[n_adds]  u32    del_p2[n_dels]  u32
 *
 * Index entries use the MSB of data_offset as a checkpoint flag
 * so readers can scan for the nearest checkpoint without
 * touching block data.
 */

namespace ibdf {

inline constexpr uint32_t MAGIC   = 0x33444249U;  // "IBD3" little-endian
inline constexpr uint16_t VERSION = 3;
inline constexpr size_t   ALIGN   = 32;

inline constexpr uint64_t CHECKPOINT_FLAG = uint64_t(1) << 63;

constexpr size_t align_up(size_t x) {
    return (x + (ALIGN - 1)) & ~(ALIGN - 1);
}

// ── On-disk structures ──────────────────────────────────────────────

struct __attribute__((packed)) FileHeader {
    uint32_t magic                = MAGIC;
    uint16_t version              = VERSION;
    uint16_t flags                = 0;
    uint64_t n_positions          = 0;
    uint64_t index_offset         = 0;
    uint32_t n_samples            = 0;
    uint32_t checkpoint_interval  = 0;
    uint8_t  _pad[32]            = {};
};
static_assert(sizeof(FileHeader) == 64);

struct __attribute__((packed)) IndexEntry {
    uint64_t bp_pos;
    uint64_t data_offset;   // MSB set = checkpoint block

    [[nodiscard]] bool is_checkpoint() const {
        return (data_offset & CHECKPOINT_FLAG) != 0;
    }
    [[nodiscard]] uint64_t offset() const {
        return data_offset & ~CHECKPOINT_FLAG;
    }
};
static_assert(sizeof(IndexEntry) == 16);

inline constexpr uint32_t BLOCK_DELTA      = 0;
inline constexpr uint32_t BLOCK_CHECKPOINT = 1;

struct __attribute__((packed)) BlockHeader {
    uint32_t block_type;    // BLOCK_DELTA or BLOCK_CHECKPOINT
    uint32_t count_a;       // checkpoint: n_pairs,  delta: n_adds
    uint32_t count_b;       // checkpoint: 0,        delta: n_dels
};

// ── Decoded views (zero-copy into mmap'd memory) ────────────────────

struct CheckpointView {
    uint64_t              bp_pos = 0;
    std::span<const float>    cm;
    std::span<const uint32_t> p1;
    std::span<const uint32_t> p2;

    [[nodiscard]] uint32_t n_pairs() const {
        return static_cast<uint32_t>(cm.size());
    }
};

struct DeltaView {
    uint64_t              bp_pos = 0;
    std::span<const float>    add_cm;
    std::span<const uint32_t> add_p1;
    std::span<const uint32_t> add_p2;
    std::span<const float>    del_cm;
    std::span<const uint32_t> del_p1;
    std::span<const uint32_t> del_p2;

    [[nodiscard]] uint32_t n_adds() const {
        return static_cast<uint32_t>(add_cm.size());
    }
    [[nodiscard]] uint32_t n_dels() const {
        return static_cast<uint32_t>(del_cm.size());
    }
};

// ── Decoders ────────────────────────────────────────────────────────

inline const BlockHeader* peek_header(const uint8_t* base) {
    return reinterpret_cast<const BlockHeader*>(base);
}

namespace detail {
    struct ArrayReader {
        const uint8_t* p;

        std::span<const float> f32(uint32_t n) {
            auto* ptr = reinterpret_cast<const float*>(p);
            p += align_up(static_cast<size_t>(n) * 4);
            return {ptr, n};
        }
        std::span<const uint32_t> u32(uint32_t n) {
            auto* ptr = reinterpret_cast<const uint32_t*>(p);
            p += align_up(static_cast<size_t>(n) * 4);
            return {ptr, n};
        }
    };
}

inline CheckpointView decode_checkpoint(const uint8_t* base) {
    const auto* hdr = reinterpret_cast<const BlockHeader*>(base);
    uint32_t n = hdr->count_a;

    detail::ArrayReader r{base + align_up(sizeof(BlockHeader))};
    CheckpointView v;
    v.cm = r.f32(n);
    v.p1 = r.u32(n);
    v.p2 = r.u32(n);
    return v;
}

inline DeltaView decode_delta(const uint8_t* base) {
    const auto* hdr = reinterpret_cast<const BlockHeader*>(base);
    uint32_t na = hdr->count_a, nd = hdr->count_b;

    detail::ArrayReader r{base + align_up(sizeof(BlockHeader))};
    DeltaView v;
    v.add_cm = r.f32(na);  v.add_p1 = r.u32(na);  v.add_p2 = r.u32(na);
    v.del_cm = r.f32(nd);  v.del_p1 = r.u32(nd);  v.del_p2 = r.u32(nd);
    return v;
}

// ── Segment for active-set tracking ─────────────────────────────────

struct Segment {
    float    cm;
    uint32_t p1;
    uint32_t p2;

    bool operator==(const Segment& o) const {
        return p1 == o.p1 && p2 == o.p2 && cm == o.cm;
    }
    bool operator<(const Segment& o) const {
        if (p1 != o.p1) return p1 < o.p1;
        if (p2 != o.p2) return p2 < o.p2;
        return cm < o.cm;
    }
};

// ── Active set operations ───────────────────────────────────────────

inline void load_checkpoint(std::vector<Segment>& active,
                            const CheckpointView& cp) {
    active.clear();
    active.reserve(cp.n_pairs());
    for (uint32_t i = 0; i < cp.n_pairs(); ++i) {
        active.push_back({cp.cm[i], cp.p1[i], cp.p2[i]});
    }
    std::sort(active.begin(), active.end());
}

inline void apply_delta(std::vector<Segment>& active,
                        const DeltaView& dv) {
    // Remove dels via binary search
    for (uint32_t i = 0; i < dv.n_dels(); ++i) {
        Segment target{dv.del_cm[i], dv.del_p1[i], dv.del_p2[i]};
        auto it = std::lower_bound(active.begin(), active.end(), target);
        if (it != active.end() && *it == target) {
            active.erase(it);
        }
    }
    // Append new segments, sort the new tail, merge into sorted prefix
    size_t old_size = active.size();
    for (uint32_t i = 0; i < dv.n_adds(); ++i) {
        active.push_back({dv.add_cm[i], dv.add_p1[i], dv.add_p2[i]});
    }
    if (dv.n_adds() > 0) {
        std::sort(active.begin() + old_size, active.end());
        std::inplace_merge(active.begin(), active.begin() + old_size, active.end());
    }
}

} // namespace ibdf
