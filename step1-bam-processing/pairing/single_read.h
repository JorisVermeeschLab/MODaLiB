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


void print(const ReadInfo& read_info, const TandemRepeat& tr) {

    uint16_t ms_start = read_info.start + tr.start;
    uint16_t ms_end = read_info.start + tr.end;

    std::string s1;
    std::string s1a;
    int i = 0;
    for (const auto& pair : read_info.bases) {
        uint16_t position = pair.first;
        s1.append(1, pair.second.c);
        if ((!pair.second.inserted) && (!pair.second.deleted)) {
            if ((ms_start <= position) && (position <= ms_end)) {
                s1a.append(1, 'R');
            } else {
                s1a.append(1, pair.second.type_as_char());
            }
        } else {
            s1a.append(1, pair.second.type_as_char());
        }
        i++;
    }
    std::cout << s1 << " " << read_info.strand << std::endl;
    std::cout << s1a << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
}


void update_statistics_with_base(const BaseInfo& base, Fragment& fragment) {

    if (!base.uninformative) {
        fragment.s_count++;
        if (base.is_cpg_site) {
            if (base.methylated) fragment.n_methylated++;
            else fragment.n_unmethylated++;
        }
        if (base.putative_snp) fragment.snp_count++;

        // If cytosine or guanine, increase GC content, and vice versa
        if (base.is_c_or_g()) fragment.gc_count++;
        else fragment.at_count++;
    }
}


Fragment compute_fragment_statistics(const ReadInfo& read_info, const AlignmentRefSequence& ref, uint8_t min_phred) {

    TandemRepeat tr = find_microsatellite(read_info.ref_bases);
    if (tr.is_valid) {
        std::cout << read_info.id << std::endl;
        std::cout << "Found repeat motif: " << tr.motif << std::endl;
        std::cout << "Number of repeats: " << tr.n_repeats << std::endl;
        print(read_info, tr);
    }

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
        update_statistics_with_base(pair.second, fragment);
    }
    return fragment;
}

#endif //MODALIB_SINGLE_READ_H
