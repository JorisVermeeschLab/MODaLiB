#ifndef _READ_INFO_H
#define _READ_INFO_H

#include <cstdint>
#include <utility>
#include <sstream>
#include <string>

#include "deserialization.h"

namespace modalib {

static constexpr size_t FRAGMENT_SIZE_IN_BYTES = 22;

struct FragmentInfo {  // TODO: should be called FragmentInfo, or Fragment

    // Genomic coordinates of the fragment start
    uint8_t start_chr_id;
    uint32_t start_pos;

    // Genomic coordinates of the fragment end
    uint8_t end_chr_id;
    uint32_t end_pos;

    // Fragment length
    uint16_t size;

    // Number of methylated CpG sites
    uint16_t n_methylated = 0;

    // Number of non-methylated CpG sites
    uint16_t n_unmethylated = 0;

    // Number of single nucleotide polymorphisms
    uint16_t snp_count = 0;

    // Number of cytosines and guanines
    uint16_t gc_count = 0;

    // Number of informative bases
    uint16_t s_count = 0;

    // Constructor
    FragmentInfo() = default;
    FragmentInfo(const FragmentInfo&) = default;
    FragmentInfo& operator=(const FragmentInfo&) = default;

    auto operator<(FragmentInfo& other) const {
    	if (start_chr_id < other.start_chr_id) return true;
    	else if (start_chr_id > other.start_chr_id) return false;
    	else return (start_pos < other.start_pos);
    }
    inline int end() const { return start_pos + size - 1; }
    inline int n_sites() const { return n_methylated + n_unmethylated; }
    inline float get_gc_ratio() const {
        if (s_count == 0) return 0.0f;
        else return static_cast<float>(gc_count) / static_cast<float>(s_count);
    }
    inline float get_meth_ratio() const {
        return static_cast<float>(n_methylated) / (static_cast<float>(n_methylated + n_unmethylated) + 0.001f);
    }
    inline float get_mutation_burden() const {
        if (s_count == 0) return 0.0f;
        else return static_cast<float>(snp_count) / static_cast<float>(s_count);
    }
    inline bool chromosomes_mismatch() const {
        return start_chr_id != end_chr_id;
    }
    void deserialize(const char* buffer) {
        int offset = 0;
        offset += Deserializer::deserialize(buffer + offset, this->start_chr_id);
        offset += Deserializer::deserialize(buffer + offset, this->start_pos);
        offset += Deserializer::deserialize(buffer + offset, this->end_chr_id);
        offset += Deserializer::deserialize(buffer + offset, this->end_pos);
        offset += Deserializer::deserialize(buffer + offset, this->size);
        offset += Deserializer::deserialize(buffer + offset, this->n_methylated);
        offset += Deserializer::deserialize(buffer + offset, this->n_unmethylated);
        offset += Deserializer::deserialize(buffer + offset, this->snp_count);
        offset += Deserializer::deserialize(buffer + offset, this->gc_count);
        offset += Deserializer::deserialize(buffer + offset, this->s_count);

        // Post-processing
        if ((start_chr_id == end_chr_id) && (start_pos > end_pos)) {
            uint32_t pos = start_pos;
            this->start_pos = end_pos;
            this->end_pos = pos;
        }
    }
    std::string to_string() const {
        std::stringstream ss;
        ss << "Fragment(";
        ss << "chr" << std::to_string(start_chr_id) << " " << std::to_string(start_pos) << ", ";
        ss << "chr" << std::to_string(end_chr_id) << " " << std::to_string(end_pos) << ", ";
        ss << std::to_string(size) << ", ";
        ss << std::to_string(n_methylated) << ", ";
        ss << std::to_string(n_unmethylated) << ", ";
        ss << std::to_string(snp_count) << ", ";
        ss << std::to_string(gc_count) << ", ";
        ss << std::to_string(s_count);
        ss << ")";
        return ss.str();
    }
};

} // namespace modalib

#endif //_READ_INFO_H
