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

#ifndef MODALIB_FRAGMENT_INFO_H
#define MODALIB_FRAGMENT_INFO_H

#include <cstdint>
#include <utility>
#include <sstream>
#include <string>
#include <cmath>
#include <assert.h>

#include "deserialization.h"

namespace modalib {

// Number of bytes used to represent each fragment in the input binary file
static constexpr size_t FRAGMENT_SIZE_IN_BYTES = 29;

struct FragmentInfo {

    // Genomic coordinates
    uint8_t chr_id = 0;
    uint32_t start_pos = 0;
    uint32_t end_pos = 0;

    // Fragment length
    uint16_t size = 0 ;
    uint16_t insert_size = 0;

    // Methylation statuses
    uint8_t n_methylated_cpg_r1 = 0;
    uint8_t n_unmethylated_cpg_r1 = 0;
    uint8_t n_methylated_cpg_r2 = 0;
    uint8_t n_unmethylated_cpg_r2 = 0;
    uint8_t n_methylated_cpg = 0;
    uint8_t n_unmethylated_cpg = 0;
    uint8_t n_methylated_ch = 0;
    uint8_t n_unmethylated_ch = 0;

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

    // End motifs
    uint8_t end_motif_5p = 0;
    uint8_t end_motif_3p = 0;

    // Mapping quality
    uint8_t mapping_quality = 0;

    // Strand
    uint8_t strand = 0;

    // Fragment weight
    float weight = 1.0f;

    // Constructor
    FragmentInfo() = default;
    FragmentInfo(const FragmentInfo&) = default;
    FragmentInfo& operator=(const FragmentInfo&) = default;

    auto operator<(FragmentInfo& other) const {
    	if (chr_id < other.chr_id) return true;
    	else if (chr_id > other.chr_id) return false;
    	else return (start_pos < other.start_pos);
    }
    inline int n_cpg_sites() const { return n_methylated_cpg + n_unmethylated_cpg; }
    inline int n_ch_sites() const { return n_methylated_ch + n_unmethylated_ch; }
    inline float get_gc_ratio() const {
        if (s_count == 0) return 0.0f;
        else return static_cast<float>(gc_count) / static_cast<float>(s_count);
    }
    inline float get_r1_meth_ratio() const {
        auto n = static_cast<float>(n_methylated_cpg_r1 + n_unmethylated_cpg_r1);
        if (n == 0) return 0.0f;
        else return static_cast<float>(n_methylated_cpg_r1) / n;
    }
    inline float get_r2_meth_ratio() const {
        auto n = static_cast<float>(n_methylated_cpg_r2 + n_unmethylated_cpg_r2);
        if (n == 0) return 0.0f;
        else return static_cast<float>(n_methylated_cpg_r2) / n;
    }
    inline float get_cpg_meth_ratio() const {
        auto n = static_cast<float>(n_methylated_cpg + n_unmethylated_cpg);
        if (n == 0) return 0.0f;
        else return static_cast<float>(n_methylated_cpg) / n;
    }
    inline float get_ch_meth_ratio() const {
        auto n = static_cast<float>(n_methylated_ch + n_unmethylated_ch);
        if (n == 0) return 0.0f;
        else return static_cast<float>(n_methylated_ch) / n;
    }
    inline float get_mutation_burden() const {
        if (s_count == 0) return 0.0f;
        else return static_cast<float>(snp_count) / static_cast<float>(s_count);
    }
    inline float overlap(const FragmentInfo& other) const {
        if (this->chr_id != other.chr_id) {
            return 0.0f;
        }
        // TODO: what about this->chr_id != other.end_chr_id (e.g., translocations)?
        float min1 = this->start_pos;
        float max1 = this->start_pos + this->size - 1;
        float min2 = other.start_pos;
        float max2 = other.start_pos + other.size - 1;
        float numerator = std::fmax(0.0f, std::fmin(max1, max2) - std::fmax(min1, min2) + 1);
        float denominator = std::fmax(max1 - min1 + 1, max2 - min2 + 1);
        if (denominator == 0.0f) return 0.0f;
        return numerator / denominator;
    }
    void serialize(char* buffer) const {
        int offset = 0;
        offset += Deserializer::serialize(buffer + offset, this->chr_id);
        offset += Deserializer::serialize(buffer + offset, this->start_pos);
        offset += Deserializer::serialize(buffer + offset, this->insert_size);
        offset += Deserializer::serialize(buffer + offset, this->size);
        offset += Deserializer::serialize(buffer + offset, this->n_methylated_cpg_r1);
        offset += Deserializer::serialize(buffer + offset, this->n_unmethylated_cpg_r1);
        offset += Deserializer::serialize(buffer + offset, this->n_methylated_cpg_r2);
        offset += Deserializer::serialize(buffer + offset, this->n_unmethylated_cpg_r2);
        offset += Deserializer::serialize(buffer + offset, this->n_methylated_cpg);
        offset += Deserializer::serialize(buffer + offset, this->n_unmethylated_cpg);
        offset += Deserializer::serialize(buffer + offset, this->n_methylated_ch);
        offset += Deserializer::serialize(buffer + offset, this->n_unmethylated_ch);
        offset += Deserializer::serialize(buffer + offset, this->gc_count);
        offset += Deserializer::serialize(buffer + offset, this->s_count);
        offset += Deserializer::serialize(buffer + offset, this->end_motif_5p);
        offset += Deserializer::serialize(buffer + offset, this->end_motif_3p);
        offset += Deserializer::serialize(buffer + offset, this->mapping_quality);
        offset += Deserializer::serialize(buffer + offset, this->strand);
        offset += Deserializer::serialize(buffer + offset, this->weight);
        assert(offset == FRAGMENT_SIZE_IN_BYTES);
    }
    void deserialize(const char* buffer) {
        int offset = 0;
        offset += Deserializer::deserialize(buffer + offset, this->chr_id);
        offset += Deserializer::deserialize(buffer + offset, this->start_pos);
        offset += Deserializer::deserialize(buffer + offset, this->insert_size);
        offset += Deserializer::deserialize(buffer + offset, this->size);
        this->end_pos = this->start_pos + this->size - 1;
        offset += Deserializer::deserialize(buffer + offset, this->n_methylated_cpg_r1);
        offset += Deserializer::deserialize(buffer + offset, this->n_unmethylated_cpg_r1);
        offset += Deserializer::deserialize(buffer + offset, this->n_methylated_cpg_r2);
        offset += Deserializer::deserialize(buffer + offset, this->n_unmethylated_cpg_r2);
        offset += Deserializer::deserialize(buffer + offset, this->n_methylated_cpg);
        offset += Deserializer::deserialize(buffer + offset, this->n_unmethylated_cpg);
        offset += Deserializer::deserialize(buffer + offset, this->n_methylated_ch);
        offset += Deserializer::deserialize(buffer + offset, this->n_unmethylated_ch);
        offset += Deserializer::deserialize(buffer + offset, this->gc_count);
        offset += Deserializer::deserialize(buffer + offset, this->s_count);
        offset += Deserializer::deserialize(buffer + offset, this->end_motif_5p);
        offset += Deserializer::deserialize(buffer + offset, this->end_motif_3p);
        offset += Deserializer::deserialize(buffer + offset, this->mapping_quality);
        offset += Deserializer::deserialize(buffer + offset, this->strand);
        offset += Deserializer::deserialize(buffer + offset, this->weight);
        assert(offset == FRAGMENT_SIZE_IN_BYTES);
    }
    std::string to_string() const {
        std::stringstream ss;
        ss << "Fragment(";
        ss << "chr" << std::to_string(chr_id) << " " << std::to_string(start_pos) << "-";
        ss << std::to_string(end_pos) << ", ";
        ss << std::to_string(size) << ", ";
        ss << std::to_string(n_methylated_cpg) << ", ";
        ss << std::to_string(n_unmethylated_cpg) << ", ";
        ss << std::to_string(snp_count) << ", ";
        ss << std::to_string(gc_count) << ", ";
        ss << std::to_string(s_count);
        ss << ")";
        return ss.str();
    }
};

} // namespace modalib

#endif  // MODALIB_FRAGMENT_INFO_H
