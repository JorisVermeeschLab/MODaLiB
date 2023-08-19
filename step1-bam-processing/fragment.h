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

#include <cstring>
#include <cstdlib>
#include <stdexcept>

#include "read_info.h"
#include "serialization.h"


static constexpr size_t FRAGMENT_SIZE_IN_BYTES = 22;


struct Fragment {
public:

    // Whether fragment is single-read
    bool single_read = false;

    // Genomic coordinates of the fragment start
    uint8_t start_chr_id;
    uint32_t start_pos;

    // Genomic coordinates of the fragment end
    uint8_t end_chr_id;
    uint32_t end_pos;

    // Fragment length
    uint16_t length;

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

    // Reads
    ReadInfo r1;
    ReadInfo r2;

public:
    Fragment(ReadInfo& r) {
        this->single_read = true;
        this->r1 = r;
        this->r2 = r;  // Single-read
        this->start_chr_id = r.chr_id;
        this->start_pos = r.start;
        this->end_chr_id = r.chr_id;
        this->end_pos = r.start + r.length - 1;
        this->length = r.length;
        this->n_methylated = r.n_methylated;
        this->n_unmethylated = r.n_unmethylated;
        this->snp_count = r.snp_count;
        this->gc_count = r.gc_count;
        this->s_count = r.s_count;
    }
    Fragment(ReadInfo& r1, ReadInfo& r2) {
        this->single_read = false;
        this->r1 = r1;
        this->r2 = r2;
        this->start_chr_id = r1.chr_id;
        this->start_pos = r1.start;
        this->end_chr_id = r2.chr_id;
        this->end_pos = r2.start + r2.length - 1;
        this->length = r1.insert_size;
        this->n_methylated = r1.n_methylated_nov + r2.n_methylated;
        this->n_unmethylated = r1.n_unmethylated_nov + r2.n_unmethylated;
        this->snp_count = r1.snp_count_nov + r2.snp_count;
        this->gc_count = r1.gc_count_nov + r2.gc_count;
        this->s_count = r1.s_count_nov + r2.s_count;
    }
    bool failed_quality_checks() {
        return this->r1.failed_quality_checks() && this->r2.failed_quality_checks();
    }
    bool is_duplicate() {
        return this->r1.is_duplicate() && this->r2.is_duplicate();
    }
    bool mapped_in_proper_pair() {
        if (this->single_read) {
            return true;
        } else {
            return this->r1.mapped_in_proper_pair() && this->r2.mapped_in_proper_pair();
        }
    }
    void serialize(char* buffer) const {
        int offset = 0;
        offset += Serializer::serialize(buffer + offset, (uint8_t) this->start_chr_id);
        offset += Serializer::serialize(buffer + offset, (uint32_t) this->start_pos);
        offset += Serializer::serialize(buffer + offset, (uint8_t) this->end_chr_id);
        offset += Serializer::serialize(buffer + offset, (uint32_t) this->end_pos);
        offset += Serializer::serialize(buffer + offset, (uint16_t) this->length);
        offset += Serializer::serialize(buffer + offset, (uint16_t) this->n_methylated);
        offset += Serializer::serialize(buffer + offset, (uint16_t) this->n_unmethylated);
        offset += Serializer::serialize(buffer + offset, (uint16_t) this->snp_count);
        offset += Serializer::serialize(buffer + offset, (uint16_t) this->gc_count);
        Serializer::serialize(buffer + offset, (uint16_t) this->s_count);
    }
};

#endif // FRAGMENT_H_
