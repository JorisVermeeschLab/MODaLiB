//  single_read.h
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

#ifndef MODALIB_SINGLE_READ_H
#define MODALIB_SINGLE_READ_H

#include <cstdint>

#include "../read_info.h"
#include "../base.h"
#include "../fragment.h"


void print(const ReadInfo& read_info, const AlignmentRefSequence& ref) {
    std::string s1;
    std::string s1a;
    for (const auto& pair : read_info.bases) {
        s1.append(1, pair.second.c);
        s1a.append(1, pair.second.type_as_char());
    }
    std::cout << s1 << " " << read_info.strand << std::endl;
    std::cout << s1a << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
}


Fragment compute_fragment_statistics(const ReadInfo& read_info, const AlignmentRefSequence& ref, uint8_t min_phred) {

    print(read_info, ref);

    Fragment fragment;
    fragment.single_read = true;
    fragment.r1 = read_info;
    fragment.r2 = read_info;  // Single-read
    fragment.start_chr_id = read_info.chr_id;
    fragment.start_pos = read_info.start;
    fragment.end_chr_id = read_info.chr_id;
    fragment.end_pos = read_info.start + read_info.length - 1;
    fragment.length = read_info.length;

    for (const auto& pair : read_info.bases) {

        // uint16_t position = pair.first;
        const BaseInfo& base = pair.second;

        if ((base.type != BASE_TYPE_DELETED) and (base.phred_score >= min_phred)) {
            fragment.s_count++;

            if (base.type == BASE_TYPE_METHYLATED_CPG) fragment.n_methylated++;
            else if (base.type == BASE_TYPE_UNMETHYLATED_CPG) fragment.n_unmethylated++;
            else if (base.type == BASE_TYPE_PUTATIVE_SNP) fragment.snp_count++;

            // If cytosine or guanine, increase GC content, and vice versa
            if (base.is_c_or_g()) fragment.gc_count++;
            else fragment.at_count++;
        }
    }
    return fragment;
}

#endif //MODALIB_SINGLE_READ_H
