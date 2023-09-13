//  mm_ml.h
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

#ifndef MODALIB_MM_ML_H
#define MODALIB_MM_ML_H

#include <sstream>
#include <string>
#include <utility>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <utility>
#include <cassert>

#include "../read_info.h"
#include "../reference.h"
#include "../strand.h"

#include "htslib/sam.h"
#include "htslib/faidx.h"


// End of sequence
static constexpr int END_OF_SEQUENCE = -2;


std::vector<std::string> split(std::string s, const std::string& delimiter) {
    std::vector<std::string> elements;
    size_t pos;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        std::string token = s.substr(0, pos);
        elements.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    if (!s.empty()) elements.push_back(s);
    return elements;
}

inline char get_base(const uint8_t* seq, int pos, int seq_len, bool reverse) {
    if ((pos < 0) or (pos >= seq_len)) return -1;
    if (reverse) pos = seq_len - 1 - pos;
    return seq_nt16_str[bam_seqi(seq, pos)];
}

inline bool is_cytosine(char base, bool complement) {
    if (!complement) return base == 'C';
    else return base == 'G';
}

int move_to_next_cytosine(uint8_t* seq, int offset, int seq_len, bool reverse_complement) {
    offset++;
    if ((seq_len < 1) or (offset >= seq_len)) return END_OF_SEQUENCE;
    while ((offset < seq_len) && (!is_cytosine(get_base(seq, offset, seq_len, reverse_complement), reverse_complement))) {
        offset++;
    }
    if (offset >= seq_len) return END_OF_SEQUENCE;
    return offset;
}

void methylation_from_aux(bam1_t* b, ReadInfo& read, uint8_t ml_cutoff) {

    // TODO: account for hard clipping (H):
    // https://github.com/samtools/hts-specs/issues/646
    // This can be done by identifying the CpG sites and using the C-CpG pattern
    // identify the Cs that are not in the clipped part of the sequence.

    int strand = get_strand(b);

    uint8_t* ml = bam_aux_get(b, "ML");
    char* mm = bam_aux2Z(bam_aux_get(b, "MM"));

    uint8_t* seq = bam_get_seq(b);
    int seq_len = b->core.l_qseq;

    std::vector<int8_t> meth_statuses;

    int seq_pos = -1;
    std::vector<std::string> modifications = split(mm, ";");

    size_t ml_pos = 0;

    for (const std::string& modification : modifications) {

        // TODO: use regex to parse "C+m", "C-m", "C+hm", etc.
        if ((modification.rfind("C+m", 0) == 0) or (modification.rfind("C-m", 0) == 0)) {

            int mm_strand = (modification.rfind("C-m", 0) == 0) ? 0 : 1;
            bool reverse_complement = (strand & 1) ^ mm_strand;

            // If '?' is present in the operation, then it means that
            // we don't know the methylation status for skipped cytosines.
            // If '.' is present instead, then the skipped cytosines
            // are indeed unmodified (unmethylated).
            std::vector<std::string> elements = split(modification, ",");
            bool skipped_are_unknown = (elements.at(0).back() == '?');

            for (int i = 1; i < elements.size(); i++) {

                // Skip unmodified cytosines, or treat them as unmethylated (if '.' mode)
                for (int k = 0; k < std::stoi(elements.at(i)); k++) {
                    seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, reverse_complement);
                    if (!skipped_are_unknown) meth_statuses.push_back(0);
                    else meth_statuses.push_back(-1);
                }

                // Go to next modified cytosine. If confidence exceeds the threshold,
                // base is considered as methylated and vice versa.
                seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, reverse_complement);
                if (ml[ml_pos] >= ml_cutoff) meth_statuses.push_back(1);
                else meth_statuses.push_back(0);

                // Go to next operation
                ml_pos++;
            }

            // Treat remaining unmodified cytosines as unmethylated (if '.' mode)
            while (seq_pos != END_OF_SEQUENCE) {
                seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, reverse_complement);
                if (seq_pos == END_OF_SEQUENCE) break;
                if (!skipped_are_unknown) meth_statuses.push_back(0);
                else meth_statuses.push_back(-1);
            }

            break;
        } else {
            std::vector<std::string> elements = split(modification, ",");
            ml_pos += (!elements.empty()) ? elements.size() - 1 : 0;
        }
    }

    // Update read statistics
    int n_cytosines = 0;
    int i = 0;
    for (auto& pair : read.bases) {
        BaseInfo& base = pair.second;
        if (is_cytosine(base.c, (strand & 1) == 0)) {
            if (i >= meth_statuses.size()) break;  // TODO: show warning message
            int8_t status = meth_statuses.at(i++);
            if (status == 1) base.methylated = 1;
            if (status == 0) base.methylated = 0;
            n_cytosines++;
        }
    }
    assert(n_cytosines == meth_statuses.size());
}

#endif // MODALIB_MM_ML_H
