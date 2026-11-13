//  fragment.h
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

#ifndef FRAGMENT_H_
#define FRAGMENT_H_

#include <algorithm>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>

#include "strand.h"
#include "read_info.h"
#include "serialization.h"
#include "msi/vntr.h"
#include "end_motifs.h"
#include "methylation/utils.h"


// Number of bytes used to represent a fragment in the output binary format
static constexpr size_t FRAGMENT_SIZE_IN_BYTES = 29;


struct Fragment {
public:

    // Whether fragment is single-read
    bool single_read = false;

    // Genomic coordinates
    uint8_t chr_id = 0;
    uint32_t start_pos = 1;
    uint32_t end_pos = 1;

    // Fragment length
    uint16_t length = 0;

    // Insert size
    uint16_t insert_size = 0;

    // List of cytosines
    std::vector<std::pair<uint32_t, BaseInfo>> cytosines;

    // Methylation of individual reads
    uint16_t n_methylated_cpg_r1 = 0;
    uint16_t n_unmethylated_cpg_r1 = 0;
    uint16_t n_methylated_cpg_r2 = 0;
    uint16_t n_unmethylated_cpg_r2 = 0;

    // Methylation in CpG contexts
    uint16_t n_methylated_cpg = 0;
    uint16_t n_unmethylated_cpg = 0;

    // Methylation in CHG contexts
    uint16_t n_methylated_chg = 0;
    uint16_t n_unmethylated_chg = 0;

    // Methylation in CHH contexts
    uint16_t n_methylated_chh = 0;
    uint16_t n_unmethylated_chh = 0;

    // Number of single nucleotide polymorphisms
    uint16_t snp_count = 0;

    // Number of cytosines and guanines
    uint16_t gc_count = 0;

    // Number of informative bases = gc_count + at_count
    uint16_t s_count = 0;

    // Number of inserted bases
    uint16_t n_inserted = 0;

    // Number of deleted bases
    uint16_t n_deleted = 0;

    // Number of soft-clipped bases
    uint16_t  n_soft_clipped = 0;

    // Strand
    int strand = STRAND_OT;

    // End motifs
    EndMotif end_motif_5p;
    EndMotif end_motif_3p;

    // Microsatellites
    std::vector<TandemRepeat> mss;

    // Reads
    ReadInfo r1;
    ReadInfo r2;

    // Weight
    float weight = 1.0f;

public:

    Fragment() = default;
    ~Fragment() = default;

    uint8_t mapping_quality() const {
        auto qual1 = static_cast<int>(r1.mapq);
        auto qual2 = static_cast<int>(r2.mapq);
        return (qual1 + qual2) / 2;
    }
    bool failed_quality_checks() const {
        return this->r1.failed_quality_checks() && this->r2.failed_quality_checks();
    }
    bool is_duplicate() const {
        return this->r1.is_duplicate() && this->r2.is_duplicate();
    }
    bool is_valid() const {
        if (!single_read) {
            if (r1.is_forward_strand() == r2.is_forward_strand()) return false;
            if (r1.strand != r2.strand) return false;
        }
        if (length == 0) return false;
        if (failed_quality_checks()) return false;
        if (is_duplicate()) return false;
        if (chr_id == 0) return false;
        if (!end_motif_5p.is_valid()) return false;
        return true;
    }
    void serialize(char* buffer) const {

        // Length of the reference sequence in the alignment
        uint16_t r_len = this->insert_size;
        // assert(this->start_pos + r_len == this->end_pos + 1);
        if (r1.chr_id == r2.chr_id) assert(this->start_pos <= this->end_pos);

        // Length of the cfDNA fragment in the alignment
        uint16_t q_len = this->length;

        // Number of methylated/unmethylated cytosines
        auto cpg_meth_counts_r1 = methylation_counts_as_uint8_t(this->n_methylated_cpg_r1, this->n_unmethylated_cpg_r1);
        auto cpg_meth_counts_r2 = methylation_counts_as_uint8_t(this->n_methylated_cpg_r2, this->n_unmethylated_cpg_r2);
        auto cpg_meth_counts = methylation_counts_as_uint8_t(this->n_methylated_cpg, this->n_unmethylated_cpg);
        auto ch_meth_counts = methylation_counts_as_uint8_t(
            this->n_methylated_chg + this->n_methylated_chh,
            this->n_unmethylated_chg + this->n_unmethylated_chh
        );

        // Serialize fragment info
        int offset = 0;
        offset += Serializer::serialize(buffer + offset, (uint8_t) this->chr_id);
        offset += Serializer::serialize(buffer + offset, (uint32_t) this->start_pos);
        offset += Serializer::serialize(buffer + offset, (uint16_t) r_len);
        offset += Serializer::serialize(buffer + offset, (uint16_t) q_len);
        offset += Serializer::serialize(buffer + offset, (uint8_t) cpg_meth_counts_r1.first);
        offset += Serializer::serialize(buffer + offset, (uint8_t) cpg_meth_counts_r1.second);
        offset += Serializer::serialize(buffer + offset, (uint8_t) cpg_meth_counts_r2.first);
        offset += Serializer::serialize(buffer + offset, (uint8_t) cpg_meth_counts_r2.second);
        offset += Serializer::serialize(buffer + offset, (uint8_t) cpg_meth_counts.first);
        offset += Serializer::serialize(buffer + offset, (uint8_t) cpg_meth_counts.second);
        offset += Serializer::serialize(buffer + offset, (uint8_t) ch_meth_counts.first);
        offset += Serializer::serialize(buffer + offset, (uint8_t) ch_meth_counts.second);
        offset += Serializer::serialize(buffer + offset, (uint16_t) this->gc_count);
        offset += Serializer::serialize(buffer + offset, (uint16_t) this->s_count);
        offset += Serializer::serialize(buffer + offset, (uint8_t) this->end_motif_5p.as_uint8());
        offset += Serializer::serialize(buffer + offset, (uint8_t) this->end_motif_3p.as_uint8());
        offset += Serializer::serialize(buffer + offset, (uint8_t) this->mapping_quality());
        offset += Serializer::serialize(buffer + offset, (uint8_t) this->strand);
        offset += Serializer::serialize(buffer + offset, (float) this->weight);

        assert(FRAGMENT_SIZE_IN_BYTES == offset);
    }
};

#endif // FRAGMENT_H_
