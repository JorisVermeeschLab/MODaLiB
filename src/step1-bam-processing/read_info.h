//  read_info.h
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

#ifndef READ_INFO_H_
#define READ_INFO_H_

#include <cstring>
#include <cstdlib>
#include <unordered_set>
#include <vector>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

#include "base.h"
#include "params.h"


/**
 * Parse chromosome name and return chromosome identifier.
 *
 * @param c Chromosome name.
 * @return Chromosome identifier.
 */
inline uint8_t parse_chr_name(const char* c) {
    size_t n = strlen(c);
    uint8_t chr;

    // If string starts with "chr", ignore the 3 first characters
    if ((n > 3) && (c[0] == 'c') && (c[1] == 'h') && (c[2] == 'r')) {
        c += 3;
        n -= 3;
    }

    char *p;
    long converted = strtol(c, &p, 10);
    if (*p) { // If not numeric
        if (n == 1) {
            if (c[0] == 'X') chr = 23;
            else if (c[0] == 'Y') chr = 24;
            else if (c[0] == 'M') chr = 25;
            else chr = 0;
        } else if ((n == 3) && (strncmp(c, "EBV", 3) == 0)) {
            chr = 26;
        } else if (n == 2) {
            if ((c[0] == 'M') && (c[1] == 'T')) chr = 25;
            else chr = 0;
        } else {
            chr = 0;
        }
    }
    else { // If numeric
        chr = (uint8_t) converted;
    }
    return chr;
}


struct ReadInfo {
public:

    // Read identifier
    std::string id;

    // Chromosome id
    uint8_t chr_id;

    // Left-most and right-most coordinates on chromosome (reference genome).
    // Flanking soft-clipped parts have been disregarded.
    uint32_t start;
    uint32_t end;

    // Insert size
    uint16_t insert_size;

    // Read strand
    int strand;

    // Read mapping quality
    uint8_t mapq;

    // Average Phred score
    uint8_t avg_phred;

    // SAM flag
    uint16_t sam_flag;

    // Weight
    float weight = 1.0f;

    // List of base-level statistics.
    // NOTE: Soft-clipped bases are not listed here.
    std::vector<std::pair<uint32_t, BaseInfo>> bases;
    std::vector<std::pair<uint32_t, BaseInfo>> ref_bases;

public:
    ReadInfo() = default;
    ReadInfo(std::string& id, const char* chr_name, uint32_t start, uint32_t end,
             uint16_t insert_size, int strand, uint8_t mapq, uint16_t sam_flag) :
        id(id),
        chr_id(parse_chr_name(chr_name)),
        start(start),
        end(end),
        insert_size(insert_size),
        strand(strand),
        mapq(mapq),
        avg_phred(255),
        sam_flag(sam_flag) {}
    ReadInfo(const ReadInfo& other) = default;
    void trim_alignment_inplace_(const Parameters& params, bool left);
    void trim_alignment_inplace(const Parameters& params);
    inline bool is_first_in_pair() const { return sam_flag & BAM_FREAD1; }
    inline bool is_forward_strand() const { return !is_reverse_strand(); }
    inline bool is_reverse_strand() const { return sam_flag & BAM_FREVERSE; }
    inline bool mapped_in_proper_pair() const { return sam_flag & BAM_FPROPER_PAIR; }
    inline bool failed_quality_checks() const { return sam_flag & BAM_FQCFAIL; }
    inline bool is_duplicate() const { return sam_flag & BAM_FDUP; }
    inline bool is_unmapped() const { return sam_flag & BAM_FUNMAP; }
    inline bool is_secondary_alignment() const { return sam_flag & BAM_FSECONDARY; }
    inline bool is_supplementary_alignment() const { return sam_flag & BAM_FSUPPLEMENTARY; }
    inline bool paired_with(const ReadInfo& other) const {
        if (this->id != other.id) return false;
        if (this->is_unmapped() || other.is_unmapped()) return false;
        if (this->insert_size != other.insert_size) return false;

        // Do as in samtools: check whether the forward strand starts at the
        // beginning of the fragment and vice versa
        if (this->is_forward_strand() && (this->start > other.end)) return false;
        if (this->is_reverse_strand() && (this->end < other.start)) return false;

        return true;
    }
    inline bool is_worth_analyzing() const {
        return !(sam_flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP));
    }
    inline void add_base(uint32_t pos, const BaseInfo& base) {
        bases.push_back(std::pair<uint32_t, BaseInfo>(pos, base));
    }
    inline void add_ref_base(uint32_t pos, const BaseInfo& base) {
        ref_bases.push_back(std::pair<uint32_t, BaseInfo>(pos, base));
    }
};

#endif // READ_INFO_H_
