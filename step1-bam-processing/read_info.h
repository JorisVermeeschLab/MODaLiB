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

#include "base.h"


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

    // Left-most coordinate on chromosome
    uint32_t start;

    // Length of the query sequence (read)
    uint16_t length;

    // Insert size
    uint16_t insert_size;

    // SAM flag
    uint16_t sam_flag;

    // Read strand
    int strand;

    // List of base-level statistics
    std::vector<std::pair<uint16_t, BaseInfo>> bases;

    // Reference sequence
    std::vector<char> ref_bases;

public:
    ReadInfo() = default;
    ReadInfo(std::string& id, const char* chr_name, uint32_t start, uint16_t length, uint16_t insert_size, int strand, uint16_t sam_flag) :
        id(id),
        chr_id(parse_chr_name(chr_name)),
        start(start),
        length(length),
        insert_size(insert_size),
        strand(strand),
        sam_flag(sam_flag) {}
    ReadInfo(const ReadInfo& other) = default;
    inline bool is_first_in_pair() const { return sam_flag & BAM_FREAD1; }
    inline bool is_reverse_strand() const { return sam_flag & BAM_FREVERSE; }
    inline bool mapped_in_proper_pair() const { return sam_flag & BAM_FPROPER_PAIR; }
    inline bool failed_quality_checks() const { return sam_flag & BAM_FQCFAIL; }
    inline bool is_duplicate() const { return sam_flag & BAM_FDUP; }
    inline bool is_unmapped() const { return sam_flag & BAM_FUNMAP; }
    inline bool is_secondary_alignment() const { return sam_flag & BAM_FSECONDARY; }
    inline bool is_supplementary_alignment() const { return sam_flag & BAM_FSUPPLEMENTARY; }
    inline bool paired_with(const ReadInfo& other) const { return this->id == other.id; }
    inline bool is_worth_analyzing() const {
        return !(sam_flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP));
    }
    inline void add_base(uint16_t pos, const BaseInfo& base) {
        bases.push_back(std::pair<uint16_t, BaseInfo>(pos, base));
    }
    inline void add_ref_base(char base) {
        ref_bases.push_back(base);
    }
};

#endif // READ_INFO_H_
