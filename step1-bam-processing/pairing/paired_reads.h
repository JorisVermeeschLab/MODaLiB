//  paired_reads.h
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

#ifndef PAIRED_READS_H
#define PAIRED_READS_H

#include <cstdint>
#include <limits>

#include "../read_info.h"
#include "../base.h"
#include "../fragment.h"
#include "single_read.h"


void print_pair(const ReadInfo& r1, const ReadInfo& r2) {

    size_t n = r1.bases.size();
    size_t m = r2.bases.size();

    size_t i = 0;
    size_t j = 0;

    std::string s1;
    std::string s1a;
    std::string s2;
    std::string s2a;

    while ((i < n) or (j < m)) {

        uint16_t min1 = (i < n) ? r1.bases[i].first : std::numeric_limits<uint16_t>::max();
        uint16_t min2 = (j < m) ? r2.bases[j].first : std::numeric_limits<uint16_t>::max();
        uint16_t k = std::min(min1, min2);
        bool use1 = (k == min1);
        bool use2 = (k == min2);

        if (use1) {
            s1.append(1, r1.bases[i].second.c);
            s1a.append(1, r1.bases[i].second.type_as_char());
        } else {
            s1.append(1, '-');
            s1a.append(1, ' ');
        }
        if (use2) {
            s2.append(1, r2.bases[j].second.c);
            s2a.append(1, r2.bases[j].second.type_as_char());
        } else {
            s2.append(1, '-');
            s2a.append(1, ' ');
        }

        if (use1) i++;
        if (use2) j++;
    }

    std::cout << s1 << " " << r1.sam_flag << " " << r1.strand << std::endl;
    std::cout << s1a << std::endl;
    std::cout << s2a << std::endl;
    std::cout << s2 << " " << r2.sam_flag << " " << r2.strand << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
}


Fragment compute_fragment_statistics(
        const ReadInfo& read1,
        const ReadInfo& read2,
        const AlignmentRefSequence& ref,
        uint8_t min_phred) {

    // Ensure r1 is actually R1, and vice versa
    ReadInfo r1, r2;
    if (read1.is_first_in_pair()) {
        r1 = read1; r2 = read2;
    } else {
        r2 = read1; r1 = read2;
    }

    TandemRepeat tr = find_microsatellite(read1.ref_bases);
    if (tr.is_valid) {
        std::cout << read1.id << std::endl;
        std::cout << "Found repeat motif: " << tr.motif << std::endl;
        std::cout << "Number of repeats: " << tr.n_repeats << std::endl;
        print_pair(r1, r2);
    }

    Fragment fragment;
    fragment.single_read = false;
    fragment.r1 = r1;
    fragment.r2 = r2;
    fragment.start_chr_id = r1.chr_id;
    fragment.start_pos = r1.start;
    fragment.end_chr_id = r2.chr_id;
    fragment.end_pos = r2.start + r2.length - 1;

    // Determine fragment length
    // TODO: account for soft-clipping?
    if (r1.insert_size == r2.insert_size) {  // Unambiguous case
        fragment.length = r1.insert_size;
    } else {  // Ambiguous case (short fragment and overlapping reads)
        fragment.length = std::min(r1.insert_size, r2.insert_size);
    }

    size_t n = r1.bases.size();
    size_t m = r2.bases.size();

    size_t i = 0;
    size_t j = 0;

    while ((i < n) or (j < m)) {

        uint16_t min1 = (i < n) ? r1.bases[i].first : std::numeric_limits<uint16_t>::max();
        uint16_t min2 = (j < m) ? r2.bases[j].first : std::numeric_limits<uint16_t>::max();
        uint16_t k = std::min(min1, min2);
        bool use1 = (k == min1);
        bool use2 = (k == min2);

        BaseInfo base;
        if (!(use1 && use2)) {  // The base is known on one strand only
            // Take the base from the only known strand
            if (use1) base = r1.bases[i].second;
            else base = r2.bases[j].second;
        } else {
            // Take the base on the strand with the highest Phred score
            if (r1.bases[i].second.phred_score >= r2.bases[j].second.phred_score) {
                base = r1.bases[i].second;
            } else {
                base = r2.bases[j].second;
            }
        }

        update_statistics_with_base(base, fragment);

        if (use1) i++;
        if (use2) j++;
    }

    std::cout << "LENGTH: " << fragment.length << " (" << fragment.start_pos << " - ";
    std::cout << fragment.end_pos << "), " << r1.id << std::endl;

    return fragment;
}

#endif //PAIRED_READS_H
