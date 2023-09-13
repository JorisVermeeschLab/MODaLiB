//  methylation.h
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

#ifndef METHYLATION_H_
#define METHYLATION_H_

#include <sstream>
#include <string>
#include <utility>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <utility>

#include "read_info.h"
#include "reference.h"
#include "strand.h"

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

inline char get_base(const uint8_t* seq, int pos, int seq_len) {
    if ((pos < 0) or (pos >= seq_len)) return -1;
    return seq_nt16_str[bam_seqi(seq, pos)];
}

inline bool is_cytosine(char base, int strand) {
    return base == 'C';  // TODO
    if (strand & 1) return base == 'C';
    else return base == 'G';
}

inline bool is_guanine(char base, int strand) {
    if (strand & 1) return base == 'G';
    else return base == 'C';
}

bool is_cpg(uint8_t* seq, int pos, int seq_len, int strand) {
    return true; // TODO
    if (pos >= seq_len - 1) return false;
    char base = get_base(seq, pos, seq_len);
    char next_base;
    if (strand & 1) next_base = get_base(seq, pos + 1, seq_len);
    else next_base = get_base(seq, pos - 1, seq_len);
    return is_cytosine(base, strand) and is_guanine(next_base, strand);
}

int move_to_next_cytosine(uint8_t* seq, int offset, int seq_len, int strand) {
    if ((seq_len < 1) or (offset >= seq_len)) return END_OF_SEQUENCE;
    offset++;
    while ((!is_cytosine(get_base(seq, offset, seq_len), strand)) and (offset < seq_len)) {
        offset++;
    }
    if (offset >= seq_len) return END_OF_SEQUENCE;
    return offset;
}

void methylation_from_aux(bam1_t* b, ReadInfo& read, uint8_t ml_cutoff) {

    // TODO: account for hard clipping (H):
    // https://github.com/samtools/hts-specs/issues/646

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

        if (modification.rfind("C+m", 0) == 0) {  // TODO: treat "C-m" case

            // If '?' is present in the operation, then it means that
            // we don't know the methylation status for skipped cytosines.
            // If '.' is present instead, then the skipped cytosines
            // are indeed unmodified (unmethylated).
            std::vector<std::string> elements = split(modification, ",");
            bool skipped_are_unknown = (elements.at(0).back() == '?');

            for (int i = 1; i < elements.size(); i++) {

                // Skip unmodified cytosines, or treat them as unmethylated (if '.' mode)
                for (int k = 0; k < std::stoi(elements.at(i)); k++) {
                    seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                    if (!skipped_are_unknown) meth_statuses.push_back(0);
                    else meth_statuses.push_back(-1);
                }

                // Go to next modified cytosine. If confidence exceeds the threshold,
                // base is considered as methylated and vice versa.
                seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                if (is_cpg(seq, seq_pos, seq_len, strand)) {
                    if (ml[ml_pos] >= ml_cutoff) meth_statuses.push_back(1);
                    else meth_statuses.push_back(0);
                }

                // Go to next operation
                ml_pos++;
            }

            // Treat remaining unmodified cytosines as unmethylated (if '.' mode)
            while (seq_pos != END_OF_SEQUENCE) {
                seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                if (!skipped_are_unknown) meth_statuses.push_back(0);
                else meth_statuses.push_back(-1);
            }

            break;
        } else {
            std::vector<std::string> elements = split(modification, ",");
            ml_pos += elements.size();
        }
    }

    // Update read statistics
    int i = 0;
    for (auto& pair : read.bases) {
        BaseInfo& base = pair.second;
        if (base.c == 'C') {
            int8_t status = meth_statuses.at(i++);
            if (status == 1) base.methylated = 1;
            if (status == 0) base.methylated = 0;
        }
    }

}

#endif // METHYLATION_H_
