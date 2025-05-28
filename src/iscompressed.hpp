//
// Created by Bohlender,Ryan James on 9/10/19.
//

#ifndef CARVAIBD_ISCOMPRESSED_HPP
#define CARVAIBD_ISCOMPRESSED_HPP

#include <fstream>

enum class CompressionType {
    gzip,
    zstd,
    uncompressed
};

/**
* @brief Checks the first two bytes of the given file path to see if they match the gzip magic bytes
* @tparam StringT Some string type, e.g., string or string_view
* @param file_path Path to a file that may be gzipped.
* @return true if the file is gzipped, otherwise false.
*/
template<class StringT>
CompressionType is_compressed(StringT file_path) {
    const uint8_t zstdref[4] = {0x28, 0xB5, 0x2F, 0xFD};
    const uint8_t gzipref[2] = {0x1F, 0x8B};
    uint8_t magic[4];
    std::ifstream isource(file_path, std::ios_base::binary);
    if (isource.good()) {
        isource.read((char *)magic, sizeof(magic));

        if(memcmp(magic, gzipref, sizeof(gzipref)) == 0) {
            // File is gzipped.
            return CompressionType::gzip;
        } else if(memcmp(magic, zstdref, sizeof(zstdref)) == 0) {
            // Magic Number: 0xFD2FB528
            // File is zstd compressed.
            return CompressionType::zstd;
        } else {
            return CompressionType::uncompressed;
        }
    } else {
        throw(std::logic_error("File doesn't exist."));
    }
}

#endif//CARVAIBD_ISCOMPRESSED_HPP
