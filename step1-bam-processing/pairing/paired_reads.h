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

        uint32_t min1 = (i < n) ? r1.bases[i].first : std::numeric_limits<uint32_t>::max();
        uint32_t min2 = (j < m) ? r2.bases[j].first : std::numeric_limits<uint32_t>::max();
        uint32_t k = std::min(min1, min2);
        bool use1 = (k == min1);
        bool use2 = (k == min2);

        if ((use1) && (r1.bases[i].second.soft_clipped)) use2 = false;
        else if ((use2) && (r2.bases[j].second.soft_clipped)) use1 = false;

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
}


Fragment compute_fragment_statistics(const ReadInfo& read1, const ReadInfo& read2) {

    // Ensure r1 is actually R1, and vice versa
    ReadInfo r1, r2;
    if (read1.is_first_in_pair()) {
        r1 = read1; r2 = read2;
    } else {
        r2 = read1; r1 = read2;
    }
    assert(r1.is_first_in_pair());
    assert(!r2.is_first_in_pair());

    /*
    TandemRepeat tr = find_microsatellite(read1.ref_bases);
    if (tr.is_valid) {
        std::cout << read1.id << std::endl;
        std::cout << "Found repeat motif: " << tr.motif << std::endl;
        std::cout << "Number of repeats: " << tr.n_repeats << std::endl;
    }
     */

    Fragment fragment;
    fragment.single_read = false;
    fragment.r1 = r1;
    fragment.r2 = r2;
    fragment.start_chr_id = r1.chr_id;
    fragment.start_pos = r1.start;
    fragment.end_chr_id = r2.chr_id;
    fragment.end_pos = r2.end;

    size_t n = r1.bases.size();
    size_t m = r2.bases.size();

    size_t i = 0;
    size_t j = 0;

    uint16_t fragment_length = 0;
    uint16_t n_inserted = 0;
    uint16_t n_deleted = 0;
    uint16_t n_soft_clipped = 0;

    while ((i < n) or (j < m)) {

        uint32_t min1 = (i < n) ? r1.bases[i].first : std::numeric_limits<uint32_t>::max();
        uint32_t min2 = (j < m) ? r2.bases[j].first : std::numeric_limits<uint32_t>::max();
        uint32_t k = std::min(min1, min2);  // Position on chromosome
        bool use1 = (k == min1);
        bool use2 = (k == min2);

        if ((use1) && (r1.bases[i].second.soft_clipped)) use2 = false;
        else if ((use2) && (r2.bases[j].second.soft_clipped)) use1 = false;

        BaseInfo b;
        if (use1 ^ use2) {  // The base is known on one strand only
            // Take the base from the only known strand
            if (use1) b = r1.bases[i].second;
            else b = r2.bases[j].second;
            if ((!b.uninformative) && (b.is_cpg_site)) {
                if (b.methylated) fragment.n_methylated++;
                else fragment.n_unmethylated++;
            }
            // Note: do not increment the GC-count here, to avoid
            // biasing GC-content upward. Base is considered as uninformative
            // as it is covered by only one read. We still call methylation though.
        } else {  // The base is known on both strands
            const BaseInfo& b1 = r1.bases[i].second;
            const BaseInfo& b2 = r2.bases[j].second;
            b = b1;

            if (b1.is_cpg_site || b2.is_cpg_site) {  // CpG site
                if ((!b2.is_cpg_site) || (b1.phred_score >= b2.phred_score)) b = b1;
                else b = b2;
                if (!b.uninformative) {
                    if (b.methylated) fragment.n_methylated++;
                    else fragment.n_unmethylated++;
                }
            } else if (b1.putative_snp || b2.putative_snp) {  // SNP or sequencing error
                if (b1.putative_snp && b2.putative_snp) {  // SNP
                    if (b1.c != b2.c) {  // Sequencing error on both R1 and R2
                        b = b1;
                        b.uninformative = 1;
                    } else if (b1.phred_score >= b2.phred_score) {
                        b = b1;
                    } else {
                        b = b2;
                    }
                    if (!b.uninformative) fragment.snp_count++;
                } else if (!b1.putative_snp) {  // Sequencing error on R2
                    b = b1;
                } else {  // Sequencing error on R1
                    b = b2;
                }
            } else if (b1.deleted || b2.deleted) {
                if (b2.is_match) b = b2;
                else b = b1;
            }

            if (!b.uninformative) {
                if (b.is_c_or_g()) fragment.gc_count++;
                fragment.s_count++;
            }
        }

        // Update fragment length estimation
        if (b.is_match || b.inserted) {
            fragment_length++;
        }
        if (b.inserted) n_inserted++;
        if (b.deleted) n_deleted++;
        if (b.soft_clipped) n_soft_clipped++;

        if (use1) i++;
        if (use2) j++;
    }

    // Compute fragment length
    if (r1.chr_id != r2.chr_id) {
        // If reads are mapped on different chromosomes, then the estimation
        // of the fragment length cannot be based on the reference coordinates.
        fragment.length = fragment_length;
    } else {
        fragment.length = r2.end - r1.start + 1;
        if (fragment.length + n_inserted >= n_deleted) {
            fragment.length += n_inserted - n_deleted;
        }
    }

    fragment.n_inserted = n_inserted;
    fragment.n_deleted = n_deleted;

    // Show results
    /*
    if (n_deleted or n_inserted or n_soft_clipped) {
        std::cout << r1.id << std::endl;
        print_pair(r1, r2);
        std::cout << "Fragment coordinates: " << fragment.start_chr_id << " " << fragment.start_pos << " ";
        std::cout << fragment.end_chr_id << " " << fragment.end_pos << std::endl;
        std::cout << "Fragment length: " << fragment.length << std::endl;
        std::cout << "Number of methylated CpG sites: " << fragment.n_methylated << std::endl;
        std::cout << "Number of unmethylated CpG sites: " << fragment.n_unmethylated << std::endl;
        std::cout << "Number of SNPs: " << fragment.snp_count << std::endl;
        std::cout << "Number of informative bases: " << fragment.s_count << std::endl;
        std::cout << "Number of insertions: " << n_inserted << std::endl;
        std::cout << "Number of deletions: " << n_deleted << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }
     */

    return fragment;
}

#endif //PAIRED_READS_H
