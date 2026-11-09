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
#include "../end_motifs.h"
#include "../jaggedness.h"
#include "../fragment.h"
#include "../msi/vntr.h"
#include "../msi/overlap.h"
#include "../msi/region_list.h"


void print(const ReadInfo& read_info, const TandemRepeat& tr) {

    uint32_t ms_start = read_info.start + tr.start;
    uint32_t ms_end = read_info.start + tr.end;

    std::string s1;
    std::string s1a;
    int i = 0;
    for (const auto& pair : read_info.bases) {
        uint32_t position = pair.first;
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


Fragment compute_fragment_statistics(const ReadInfo& read_info) {

    // Reference sequence
    std::vector<BaseInfo> ref_seq;

    Fragment fragment;
    fragment.single_read = true;
    fragment.r1 = read_info;
    fragment.r2 = read_info;  // Single-read
    fragment.chr_id = read_info.chr_id;
    fragment.start_pos = read_info.start;
    fragment.end_pos = read_info.end;
    fragment.length = read_info.insert_size;
    fragment.insert_size = read_info.insert_size;
    fragment.strand = read_info.strand;

    std::vector<uint32_t> insertions;
    std::vector<uint32_t> deletions;
    fragment.n_soft_clipped = 0;

    for (size_t i = 0; i < read_info.bases.size(); i++) {
        const auto& pair = read_info.bases[i];
        uint32_t k = pair.first;
        const BaseInfo& base = pair.second;

        // Store reference sequence
        BaseInfo ref_base = read_info.ref_bases[i].second;
        if (!(ref_base.soft_clipped || ref_base.deleted)) {
            ref_seq.push_back(ref_base);
        }

        // Check if base is soft-clipped
        if (ref_base.soft_clipped || base.soft_clipped) fragment.n_soft_clipped++;

        // Increment meth counts for jaggedness estimation (no trimming)
        if (base.is_cpg_site()) {
            increment_methylation_count(i, fragment.length, base.methylated);
        }

        if (!base.uninformative) {
            if (base.is_cpg_site()) {
                if (base.methylated) fragment.n_methylated_cpg++;
                else fragment.n_unmethylated_cpg++;
                fragment.cytosines.emplace_back(k, base);
            } else if (base.is_chg_site()) {
                if (base.methylated) fragment.n_methylated_chg++;
                else fragment.n_unmethylated_chg++;
                fragment.cytosines.emplace_back(k, base);
            } else if (base.is_chh_site()) {
                if (base.methylated) fragment.n_methylated_chh++;
                else fragment.n_unmethylated_chh++;
                fragment.cytosines.emplace_back(k, base);
            }
            if (base.putative_snp) fragment.snp_count++;

            // If cytosine or guanine, increase GC content, and vice versa
            if (base.is_c_or_g()) fragment.gc_count++;
            fragment.s_count++;
        }

        if (base.inserted) {
            fragment.n_inserted++;
            insertions.push_back(k);
        }
        if (base.deleted) {
            fragment.n_deleted++;
            deletions.push_back(k);
        }
    }

    // Store end motifs
    if (ref_seq.size() >= 4) {
        size_t p = read_info.ref_bases.size();
        if ((fragment.strand & 1) == 1) {
            fragment.end_motif_5p = EndMotif(ref_seq[0].c, ref_seq[1].c, ref_seq[2].c, ref_seq[3].c, STRAND_OT);
            fragment.end_motif_3p = EndMotif(ref_seq[p-1].c, ref_seq[p-2].c, ref_seq[p-3].c, ref_seq[p-4].c, STRAND_OT);
        } else {
            fragment.end_motif_5p = EndMotif(ref_seq[p-1].c, ref_seq[p-2].c, ref_seq[p-3].c, ref_seq[p-4].c, STRAND_OB);
            fragment.end_motif_3p = EndMotif(ref_seq[0].c, ref_seq[1].c, ref_seq[2].c, ref_seq[3].c, STRAND_OB);
        }
    } else {
        fragment.end_motif_5p = EndMotif('N', 'N', 'N', 'N', STRAND_OT);
        fragment.end_motif_3p = EndMotif('N', 'N', 'N', 'N', STRAND_OT);
    }

    // Fragment weight
    fragment.weight = read_info.weight;

    return fragment;
}

#endif //MODALIB_SINGLE_READ_H
