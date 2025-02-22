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
#include <vector>
#include <algorithm>

#include "../read_info.h"
#include "../base.h"
#include "../fragment.h"
#include "../end_motifs.h"
#include "../msi/vntr.h"
#include "single_read.h"
#include "../msi/region_list.h"


void print_pair(const ReadInfo& r1, const ReadInfo& r2);


Fragment compute_fragment_statistics(const ReadInfo& read1, const ReadInfo& read2) {

    // Aligned reference genome
    // NOTE: soft-clipped bases are excluded. Therefore, the length of this alignment
    // can be shorter than the alignment of one of the 2 individual reads.
    std::vector<BaseInfo> ref_seq;

    std::vector<BaseInfo> indices;

    // Ensure first read is mapped on forward strand
    ReadInfo r1, r2;
    if (read1.is_forward_strand()) {
        r1 = read1; r2 = read2;
    } else {
        r2 = read1; r1 = read2;
    }

    // Initialize fragment information
    Fragment fragment;
    fragment.single_read = false;  // Paired reads
    fragment.r1 = r1;
    fragment.r2 = r2;
    fragment.chr_id = r1.chr_id;
    fragment.start_pos = r1.start;
    fragment.end_pos = r2.end;
    fragment.insert_size = std::min(r1.insert_size, r2.insert_size);
    fragment.length = fragment.insert_size;
    fragment.strand = r1.strand;

    size_t n = r1.bases.size();
    size_t m = r2.bases.size();

    size_t i = 0;
    size_t j = 0;

    std::vector<uint32_t> ins1_not2;  // Number of bases that are insertions in R1 but are absent from R2
    std::vector<uint32_t> ins2_not1;  // Number of bases that are insertions in R2 but are absent from R1
    std::vector<uint32_t> del1_not2;  // Number of bases that are deletions in R1 but are absent from R2
    std::vector<uint32_t> del2_not1;  // Number of bases that are deletions in R2 but are absent from R1
    std::vector<uint32_t> ins1_mat2;  // Number of bases that are insertions in R1 but not in R2
    std::vector<uint32_t> ins2_mat1;  // Number of bases that are insertions in R2 but not in R1
    std::vector<uint32_t> del1_mat2;  // Number of bases that are deletions in R1 but not in R2
    std::vector<uint32_t> del2_mat1;  // Number of bases that are deletions in R2 but not in R1
    std::vector<uint32_t> ins1_ins2;  // Number of bases that are insertions in both reads
    std::vector<uint32_t> del1_del2;  // Number of bases that are deletions in both reads
    fragment.n_soft_clipped = 0;  // Number of soft-clipped bases

    while ((i < n) or (j < m)) {

        uint32_t min1 = (i < n) ? r1.bases[i].first : std::numeric_limits<uint32_t>::max();
        uint32_t min2 = (j < m) ? r2.bases[j].first : std::numeric_limits<uint32_t>::max();
        uint32_t k = std::min(min1, min2);  // Position on chromosome
        bool use1 = (k == min1);
        bool use2 = (k == min2);

        // Determine reference base
        BaseInfo ref_base;
        if (use1) ref_base = r1.ref_bases[i].second;
        else ref_base = r2.ref_bases[j].second;
        if (!(ref_base.soft_clipped || ref_base.deleted)) {
            ref_seq.emplace_back(ref_base);
        }

        BaseInfo b;  // Consensus base
        if (use1 ^ use2) {  // The base is known on one strand only
            // Take the base from the only known strand
            if (use1) b = r1.bases[i].second;
            else b = r2.bases[j].second;

            // Count read-level insertions and deletions
            if (b.inserted) {
                if (use1) ins1_not2.push_back(k);
                else ins2_not1.push_back(k);
                // b.uninformative = 1;  // TODO
            } else if (b.deleted) {
                if (use1) del1_not2.push_back(k);
                else del2_not1.push_back(k);
                // b.uninformative = 1;  // TODO
            }

            // Methylation calling
            if (!b.uninformative) {
                if (b.is_cpg_site()) {
                    if (b.methylated) fragment.n_methylated_cpg++;
                    else fragment.n_unmethylated_cpg++;
                    fragment.cytosines.emplace_back(k, b);
                } else if (b.is_chg_site()) {
                    if (b.methylated) fragment.n_methylated_chg++;
                    else fragment.n_unmethylated_chg++;
                    fragment.cytosines.emplace_back(k, b);
                } else if (b.is_chh_site()) {
                    if (b.methylated) fragment.n_methylated_chh++;
                    else fragment.n_unmethylated_chh++;
                    fragment.cytosines.emplace_back(k, b);
                }
            }

            // No variant calling in the absence of the second read,
            // as it would result in a tremendous amount of false positives.
            if (b.putative_snp) b.uninformative = 1;

            // Note: do not increment the GC-count here, to avoid
            // biasing GC-content upward. Base is considered as uninformative
            // as it is covered by only one read. We still call methylation though.

        } else {  // The base is known on both strands

            // assert(r1.ref_bases[i].second.c == r2.ref_bases[j].second.c);

            const BaseInfo& b1 = r1.bases[i].second;
            const BaseInfo& b2 = r2.bases[j].second;
            b = b1;

            if (b1.has_cytosine_context() || b2.has_cytosine_context()) {  // CpG|CHG|CHH context
                if (!b2.has_cytosine_context()) b = b1;
                else if (!b1.has_cytosine_context()) b = b2;
                else if (b1.phred_score >= b2.phred_score) b = b1;
                else b = b2;
                if (!b.uninformative) {
                    if (b.is_cpg_site()) {  // CpG context
                        if (b.methylated) fragment.n_methylated_cpg++;
                        else fragment.n_unmethylated_cpg++;
                        fragment.cytosines.emplace_back(k, b);
                    } else if (b.is_chg_site()) {  // CHG context
                        if (b.methylated) fragment.n_methylated_chg++;
                        else fragment.n_unmethylated_chg++;
                        fragment.cytosines.emplace_back(k, b);
                    } else if (b.is_chh_site()) {  // CHH context
                        if (b.methylated) fragment.n_methylated_chh++;
                        else fragment.n_unmethylated_chh++;
                        fragment.cytosines.emplace_back(k, b);
                    }
                }
            } else if (b1.putative_snp || b2.putative_snp) {  // SNP or sequencing error
                if (b1.putative_snp && b2.putative_snp) {  // Putative SNP
                    if (b1.c != b2.c) {  // Sequencing error on both reads
                        b = b1;
                        b.uninformative = 1;
                    } else {  // Same mismatch on both reads -> putative mutation
                        if (b2.uninformative) b = b1;
                        else if (b1.uninformative) b = b2;
                        if ((b1.phred_score >= b2.phred_score) && (!b1.uninformative)) {
                            b = b1;
                        } else if ((b2.phred_score >= b1.phred_score) && (!b2.uninformative)) {
                            b = b2;
                        }
                        if (!b.uninformative) fragment.snp_count++;
                    }
                } else if (!b1.putative_snp) {  // Sequencing error on second read
                    b = b1;
                    b.uninformative = 0;
                } else {  // Sequencing error on first read
                    b = b2;
                    b.uninformative = 0;
                }
            } else if (b1.deleted || b2.deleted) {
                // Systematically use R1, because if the deletion does not occur on the same
                // base for both reads due to alignment issues (e.g., within a microsatellite),
                // then the deletion in R1 will be followed by a non-deletion in R1, and vice versa.
                // By doing so, the non-deleted base will count only once toward the GC count and the
                // number of informative bases.
                // We could have chosen R2 instead, but the later is more prone to sequencing errors.
                b = b1;
            } else if (b1.inserted || b2.inserted) {
                b = b1;  // Same remark as for deletions
            }

            // Count read-level insertions and deletions
            if (b1.inserted && b2.inserted) {
                ins1_ins2.push_back(k);
            } else if (b1.deleted && b2.deleted) {
                del1_del2.push_back(k);
            } else if (b1.inserted && (!b2.inserted)) {
                ins1_mat2.push_back(k);
            } else if (b1.deleted && (!b2.deleted)) {
                del1_mat2.push_back(k);
            } else if (b2.inserted && (!b1.inserted)) {
                ins2_mat1.push_back(k);
            } else if (b2.deleted && (!b1.deleted)) {
                del2_mat1.push_back(k);
            } else if ((!b1.inserted) && (!b1.deleted)) {
                if (b2.inserted) ins2_mat1.push_back(k);
                else if (b2.deleted) del2_mat1.push_back(k);
            } else if ((!b2.inserted) && (!b2.deleted)) {
                if (b1.inserted) ins1_mat2.push_back(k);
                else if (b1.deleted) del1_mat2.push_back(k);
            }
        }

        // Update GC-content
        if (!b.uninformative) {
            if (b.is_c_or_g()) fragment.gc_count++;
            fragment.s_count++;
        }

        // Update number of soft-clipped bases
        if (b.soft_clipped || ref_base.soft_clipped) fragment.n_soft_clipped++;

        // Go to next base
        if (use1) i++;
        if (use2) j++;
    }

    // Fragment-level indel calling
    fragment.n_inserted = ins1_ins2.size();
    fragment.n_deleted = del1_del2.size();
    fragment.n_inserted += std::min(ins2_mat1.size(), ins1_mat2.size());
    fragment.n_deleted += std::min(del2_mat1.size(), del1_mat2.size());
    fragment.n_inserted += ins1_not2.size() + ins2_not1.size();
    fragment.n_deleted += del1_not2.size() + del2_not1.size();

    // Compute fragment length
    if (r1.chr_id != r2.chr_id) {
        // If reads are mapped on different chromosomes, then the estimation
        // of the fragment length cannot be based on the reference coordinates.
        fragment.length = fragment.insert_size;  // Trust insert size
    } else {
        fragment.length = fragment.insert_size;
        fragment.length += fragment.n_inserted - fragment.n_deleted;
    }
    if ((fragment.length == 0) && (fragment.end_pos > fragment.start_pos)) {
        fragment.length = fragment.end_pos - fragment.start_pos + 1;
    }

    // Store end motifs
    if ((r1.ref_bases.size() >= 4) && (r2.ref_bases.size() >= 4)) {
        size_t p = ref_seq.size();
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

    // Get methylation of R1 and R2
    fragment.n_methylated_cpg_r1 = 0;
    fragment.n_unmethylated_cpg_r1 = 0;
    for (auto &pair : r1.bases) {
        BaseInfo &b = pair.second;
        if (b.is_cpg_site()) {
            if (b.methylated) fragment.n_methylated_cpg_r1++;
            else fragment.n_unmethylated_cpg_r1++;
        }
    }
    fragment.n_methylated_cpg_r2 = 0;
    fragment.n_unmethylated_cpg_r2 = 0;
    for (auto &pair : r2.bases) {
        BaseInfo &b = pair.second;
        if (b.is_cpg_site()) {
            if (b.methylated) fragment.n_methylated_cpg_r2++;
            else fragment.n_unmethylated_cpg_r2++;
        }
    }
    if (!r1.is_first_in_pair()) {
        std::swap(fragment.n_methylated_cpg_r1, fragment.n_methylated_cpg_r2);
        std::swap(fragment.n_unmethylated_cpg_r1, fragment.n_unmethylated_cpg_r2);
    }

    // Fragment weight
    fragment.weight = 0.5f * (r1.weight + r2.weight);

    // Show results
    if (false) {
        std::cout << r1.id << std::endl;
        print_pair(r1, r2);
        std::cout << "Fragment coordinates: " << 1 * fragment.chr_id << " " << fragment.start_pos << " ";
        std::cout << fragment.end_pos << std::endl;
        // std::cout << "Fragment coordinates (alt): " << fragment.r2.start << " " << fragment.r1.end << std::endl;
        std::cout << "Insert size: " << fragment.insert_size << std::endl;
        std::cout << "Fragment length: " << fragment.length << std::endl;
        std::cout << "Number of methylated CpG sites: " << fragment.n_methylated_cpg << std::endl;
        std::cout << "Number of unmethylated CpG sites: " << fragment.n_unmethylated_cpg << std::endl;
        std::cout << "Number of methylated CHG sites: " << fragment.n_methylated_chg << std::endl;
        std::cout << "Number of unmethylated CHG sites: " << fragment.n_unmethylated_chg << std::endl;
        std::cout << "Number of methylated CHH sites: " << fragment.n_methylated_chh << std::endl;
        std::cout << "Number of unmethylated CHH sites: " << fragment.n_unmethylated_chh << std::endl;
        std::cout << "Number of SNPs: " << fragment.snp_count << std::endl;
        std::cout << "GC count: " << fragment.gc_count << std::endl;
        std::cout << "Number of informative bases: " << fragment.s_count << std::endl;
        std::cout << "Number of insertions: " << fragment.n_inserted << std::endl;
        std::cout << "Number of deletions: " << fragment.n_deleted << std::endl;
        std::cout << "Phred: " << 1 * fragment.r1.avg_phred << " " << 1 * fragment.r2.avg_phred << std::endl;
        std::cout << "MAPQ: " << 1 * fragment.r1.mapq << " " << 1 * fragment.r2.mapq << std::endl;
        std::cout << "Strand: " << fragment.r1.strand << " " << fragment.r2.strand << std::endl;
        std::cout << "End motifs: " << fragment.end_motif_5p.as_str() << " - " << fragment.end_motif_3p.as_str() << std::endl;
        std::cout << "Valid pair: " << 1 * fragment.is_valid() << std::endl;
        std::cout << "INDELS " << ins1_ins2.size() << " " << del1_del2.size() << " " << ins1_not2.size() << " " << del1_not2.size() << " " << ins2_not1.size() << " " << del2_not1.size() << " " << ins1_mat2.size() << " " << del1_mat2.size() << " " << ins2_mat1.size() << " " << del2_mat1.size() << std::endl;
        for (uint32_t pos : ins1_ins2) std::cout << "INS " << pos << std::endl;
        for (uint32_t pos : ins1_mat2) std::cout << "INS " << pos << std::endl;
        for (uint32_t pos : ins1_not2) std::cout << "INS " << pos << std::endl;
        for (uint32_t pos : ins2_mat1) std::cout << "INS " << pos << std::endl;
        for (uint32_t pos : ins2_not1) std::cout << "INS " << pos << std::endl;
        for (uint32_t pos : del1_del2) std::cout << "DEL " << pos << std::endl;
        for (uint32_t pos : del1_mat2) std::cout << "DEL " << pos << std::endl;
        for (uint32_t pos : del1_not2) std::cout << "DEL " << pos << std::endl;
        for (uint32_t pos : del2_mat1) std::cout << "DEL " << pos << std::endl;
        for (uint32_t pos : del2_not1) std::cout << "DEL " << pos << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }

    return fragment;
}


void print_pair(const ReadInfo& r1, const ReadInfo& r2) {

    size_t n = r1.bases.size();
    size_t m = r2.bases.size();

    size_t i = 0;
    size_t j = 0;

    std::string s_ref;
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

        // Manage soft-clipping
        if ((use1) && (r1.bases[i].second.soft_clipped)) use2 = false;
        else if ((use2) && (r2.bases[j].second.soft_clipped)) use1 = false;
        if ((!use1) && (!use2)) {  // This should never happen
            use1 = true;
            use2 = true;
        }

        // Show reference sequence
        BaseInfo ref_base;
        if (use1) ref_base = r1.ref_bases[i].second;
        else ref_base = r2.ref_bases[j].second;
        if (ref_base.deleted) s_ref.append(1, '-');
        else s_ref.append(1, ref_base.c);

        // Show read 1
        if (use1) {
            s1.append(1, r1.bases[i].second.c);
            s1a.append(1, r1.bases[i].second.type_as_char());
        } else {
            s1.append(1, '-');
            s1a.append(1, ' ');
        }

        // Show read 2
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

    // Visualize everything
    std::cout << s_ref << " " << "REF" << std::endl;
    std::cout << s1 << " " << r1.sam_flag << " " << r1.strand << std::endl;
    std::cout << s1a << std::endl;
    std::cout << s2a << std::endl;
    std::cout << s2 << " " << r2.sam_flag << " " << r2.strand << std::endl;
}


#endif //PAIRED_READS_H
