#pragma once

#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include "ibd_format.hpp"
#include "mapped_file.hpp"

class IbdfReader {
public:
    explicit IbdfReader(const char* path) : file_(path) {
        header_ = reinterpret_cast<const ibdf::FileHeader*>(file_.data());
        if (header_->magic != ibdf::MAGIC) {
            throw std::runtime_error("Bad IBDF magic (expected v3)");
        }
        if (header_->version != ibdf::VERSION) {
            throw std::runtime_error("Version mismatch: expected v3");
        }
        auto* idx_base = reinterpret_cast<const ibdf::IndexEntry*>(
            file_.data() + header_->index_offset);
        index_ = std::span<const ibdf::IndexEntry>(idx_base, header_->n_positions);
    }

    [[nodiscard]] uint64_t n_positions()         const { return header_->n_positions; }
    [[nodiscard]] uint32_t n_samples()           const { return header_->n_samples; }
    [[nodiscard]] uint32_t checkpoint_interval() const { return header_->checkpoint_interval; }

    [[nodiscard]] std::span<const ibdf::IndexEntry> index() const { return index_; }

    [[nodiscard]] const uint8_t* block_ptr(uint64_t pos_idx) const {
        return file_.data() + index_[pos_idx].offset();
    }

    [[nodiscard]] uint64_t find_checkpoint_before(uint64_t pos_idx) const {
        int64_t i = static_cast<int64_t>(pos_idx);
        while (i >= 0 && !index_[i].is_checkpoint()) --i;
        if (i < 0) throw std::runtime_error("No checkpoint found before index");
        return static_cast<uint64_t>(i);
    }

    [[nodiscard]] std::vector<ibdf::Segment>
    reconstruct_active_set(uint64_t target_idx) const {
        uint64_t ckpt_idx = find_checkpoint_before(target_idx);

        std::vector<ibdf::Segment> active;
        auto cp = ibdf::decode_checkpoint(block_ptr(ckpt_idx));
        cp.bp_pos = index_[ckpt_idx].bp_pos;
        ibdf::load_checkpoint(active, cp);

        for (uint64_t i = ckpt_idx + 1; i <= target_idx; ++i) {
            const auto* hdr = ibdf::peek_header(block_ptr(i));
            if (hdr->block_type == ibdf::BLOCK_CHECKPOINT) {
                auto cp2 = ibdf::decode_checkpoint(block_ptr(i));
                ibdf::load_checkpoint(active, cp2);
            } else {
                auto dv = ibdf::decode_delta(block_ptr(i));
                ibdf::apply_delta(active, dv);
            }
        }
        return active;
    }

private:
    MappedFile                          file_;
    const ibdf::FileHeader*             header_ = nullptr;
    std::span<const ibdf::IndexEntry>   index_;
};
