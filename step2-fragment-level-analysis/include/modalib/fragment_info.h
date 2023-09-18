//  fragment_info.h
//
//  Copyright 2023 Antoine Passemiers <antoine.passemiers@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//  MA 02110-1301, USA.

#ifndef MODALIB_READ_INFO_H
#define MODALIB_READ_INFO_H

#include <cstdint>
#include <utility>
#include <sstream>
#include <string>
#include <cmath>

#include "deserialization.h"

namespace modalib {

// Number of bytes used to represent each fragment in the input binary file
static constexpr size_t FRAGMENT_SIZE_IN_BYTES = 22;

struct FragmentInfo {

    // Genomic coordinates of the fragment start
    uint8_t start_chr_id = 0;
    uint32_t start_pos = 0;

    // Genomic coordinates of the fragment end
    uint8_t end_chr_id = 0;
    uint32_t end_pos = 0;

    // Fragment length
    uint16_t size = 0 ;

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

    // Number of inserted bases
    uint16_t n_inserted = 0;

    // Number of deleted bases
    uint16_t n_deleted = 0;

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
    inline float overlap(const FragmentInfo& other) const {
        if ((this->start_chr_id != other.start_chr_id) || (this->end_chr_id != other.end_chr_id)) {
            return 0.0f;
        }
        // TODO: what about this->start_chr_id != other.end_chr_id (e.g., translocations)?
        float min1 = this->start_pos;
        float max1 = this->start_pos + this->size + 1;
        float min2 = other.start_pos;
        float max2 = other.start_pos + other.size + 1;
        float numerator = std::fmax(0.0f, std::fmin(max1, max2) - std::fmax(min1, min2));
        float denominator = std::fmax(max1 - min1, max2 - min2);
        if (denominator == 0.0f) return 0.0f;
        return numerator / denominator;
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
        offset += Deserializer::deserialize(buffer + offset, this->n_inserted);
        Deserializer::deserialize(buffer + offset, this->n_deleted);

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

#endif  // MODALIB_READ_INFO_H
