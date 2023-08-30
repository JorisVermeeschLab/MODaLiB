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

    // Number of single nucleotide polymorphisms
    uint16_t snp_count = 0;
    uint16_t snp_count_nov = 0;  // Non-overlapping

    // Number of cytosines and guanines
    uint16_t gc_count = 0;
    uint16_t gc_count_nov = 0;  // Non-overlapping

    // Number of informative bases
    uint16_t s_count = 0;
    uint16_t s_count_nov = 0;  // Non-overlapping

    // Length of the query sequence (read)
    uint16_t length;

    // Insert size
    uint16_t insert_size;

    // SAM flag
    uint16_t sam_flag;

    // Number of inserted bases
    uint16_t inserted = 0;

    // Number of deleted bases
    uint16_t deleted = 0;

    // Positions of methylated CpG sites on chromosome
    std::unordered_set<uint32_t> methylated_cpg_sites;

    // Positions of unmethylated CpG sites on chromosome
    std::unordered_set<uint32_t> unmethylated_cpg_sites;

    // Positions of adenines and thymines on chromosome
    std::unordered_set<uint32_t> at_positions;

    // Positions of cytosines and guanines on chromosome
    std::unordered_set<uint32_t> gc_positions;

    // Positions of "mutations" on chromosome
    std::unordered_set<uint32_t> mut_positions;

    // Positions of informative bases on chromosome
    std::unordered_set<uint32_t> informative_positions;

public:
    ReadInfo() = default;
    ReadInfo(std::string& id, const char* chr_name, uint32_t start, uint16_t length, uint16_t insert_size, uint16_t sam_flag) :
        id(id),
        chr_id(parse_chr_name(chr_name)),
        start(start),
        length(length),
        insert_size(insert_size),
        sam_flag(sam_flag),
        methylated_cpg_sites(),
        unmethylated_cpg_sites() {}
    ReadInfo(const ReadInfo& other) = default;
    inline bool is_first_in_pair() const { return 0x40 & sam_flag; }
    inline bool is_reverse_strand() const { return 0x10 & sam_flag; }
    inline bool mapped_in_proper_pair() const { return 0x3 & sam_flag; }
    inline bool failed_quality_checks() const { return 0x200 & sam_flag; }
    inline bool is_duplicate() const { return 0x400 & sam_flag; }
    inline bool paired_with(const ReadInfo& other) const { return this->id == other.id; }
    inline bool is_non_overlapping(uint32_t read_position) const {
        if (!is_first_in_pair()) return true;
        else {
            // TODO: We assume that the 2 paired reads have the same length
            // TODO: this simplifies the computation of the overlap
            return (int) read_position < (int) insert_size - (int) length;
        }
    }
    inline void add_methylated_cpg(uint32_t pos) { methylated_cpg_sites.insert(pos); }
    inline void add_unmethylated_cpg(uint32_t pos) { unmethylated_cpg_sites.insert(pos); }
    inline int num_methylated_cpg_sites() const { return methylated_cpg_sites.size(); }
    inline int num_unmethylated_cpg_sites() const { return unmethylated_cpg_sites.size(); }
    inline void add_at_base(uint32_t pos) { at_positions.insert(pos); }
    inline void add_gc_base(uint32_t pos) { gc_positions.insert(pos); }
    inline void add_mutation(uint32_t pos) { mut_positions.insert(pos); }
    inline void add_informative_base(uint32_t pos) { informative_positions.insert(pos); }
};

#endif // READ_INFO_H_
