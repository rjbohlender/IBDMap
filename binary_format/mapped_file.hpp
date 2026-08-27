#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

class MappedFile {
public:
    explicit MappedFile(const char* path) {
        int fd = ::open(path, O_RDONLY);
        if (fd < 0) throw std::runtime_error(std::string("open: ") + path);
        struct stat st{};
        ::fstat(fd, &st);
        size_ = static_cast<size_t>(st.st_size);
        base_ = static_cast<const uint8_t*>(
            ::mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd, 0));
        ::close(fd);
        if (base_ == MAP_FAILED) throw std::runtime_error("mmap failed");
    }
    ~MappedFile() {
        if (base_ && base_ != MAP_FAILED) ::munmap(const_cast<uint8_t*>(base_), size_);
    }
    MappedFile(const MappedFile&) = delete;
    MappedFile& operator=(const MappedFile&) = delete;

    [[nodiscard]] const uint8_t* data() const { return base_; }
    [[nodiscard]] size_t size() const { return size_; }

    [[nodiscard]] bool discard(size_t offset, size_t length) const {
        if (offset > size_ || length > size_ - offset) return false;

        const long raw_page_size = ::sysconf(_SC_PAGESIZE);
        if (raw_page_size <= 0) return false;
        const size_t page_size = static_cast<size_t>(raw_page_size);

        const size_t offset_remainder = offset % page_size;
        const size_t aligned_start = offset_remainder == 0
                ? offset
                : offset + (page_size - offset_remainder);
        const size_t end = offset + length;
        const size_t aligned_end = end - (end % page_size);
        if (aligned_start >= aligned_end) return true;

        return ::madvise(
                       const_cast<uint8_t*>(base_) + aligned_start,
                       aligned_end - aligned_start,
                       MADV_DONTNEED) == 0;
    }

private:
    const uint8_t* base_ = nullptr;
    size_t size_ = 0;
};
