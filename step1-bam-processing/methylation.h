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
#include "utils.h"

#include "htslib/sam.h"
#include "htslib/faidx.h"


// End of sequence
static constexpr int END_OF_SEQUENCE = -2;


std::vector<std::string> split(std::string s, const std::string delimiter) {
    std::vector<std::string> elements;
    int pos;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        std::string token = s.substr(0, pos);
        elements.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    if (s.size() > 0) elements.push_back(s);
    return elements;
}

inline char get_base(uint8_t* seq, int pos, int seq_len, int strand) {
    if (!(strand & 1)) {
        pos = seq_len - 1 - pos;
    }
    if ((pos < 0) or (pos >= seq_len)) return -1;
    return seq_nt16_str[bam_seqi(seq, pos)];
}

inline bool is_cytosine(char base, int strand) {
    if (strand & 1) return base == 'C';
    else return base == 'G';
}

inline bool is_guanine(char base, int strand) {
    if (strand & 1) return base == 'G';
    else return base == 'C';
}

bool is_cpg(uint8_t* seq, int pos, int seq_len, int strand) {
    if (pos >= seq_len - 1) return false;
    char base = get_base(seq, pos, seq_len, strand);
    char next_base = get_base(seq, pos + 1, seq_len, strand);
    return is_cytosine(base, strand) and is_guanine(next_base, strand);
}

int move_to_next_cytosine(uint8_t* seq, int offset, int seq_len, int strand) {
    if ((seq_len < 1) or (offset >= seq_len)) return END_OF_SEQUENCE;
    offset++;
    while ((!is_cytosine(get_base(seq, offset, seq_len, strand), strand)) and (offset < seq_len)) {
        offset++;
    }
    if (offset >= seq_len) return END_OF_SEQUENCE;
    return offset;
}

void methylation_from_aux(bam1_t *b, ReadInfo& read, uint8_t ml_cutoff) {

    // TODO: account for hard clipping (H):
    // https://github.com/samtools/hts-specs/issues/646

    int strand = get_strand(b);

    uint8_t* ml = bam_aux_get(b, "ML");
    char* mm = bam_aux2Z(bam_aux_get(b, "MM"));

    uint8_t* seq = bam_get_seq(b);
    int seq_len = b->core.l_qseq;

    int seq_pos = -1;
    std::vector<std::string> modifications = split(mm, ";");

    int ml_pos = 0;

    for (std::string modification : modifications) {

        if (modification.rfind("C+m", 0) == 0) {

            std::vector<std::string> elements = split(modification, ",");

            for (int i = 1; i < elements.size(); i++) {

                // Skip unmodified cytosines
                for (int k = 0; k < std::stoi(elements.at(i)); k++) {
                    seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                    // Do not call methylation as these cytosines are inconclusive
                }

                // Go to next modified cytosine. If confidence exceeds the threshold,
                // base is considered as methylated and vice versa.
                seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                if (is_cpg(seq, seq_pos, seq_len, strand)) {
                    if (ml[ml_pos] >= ml_cutoff) read.add_methylated_cpg(seq_pos);
                    else read.add_unmethylated_cpg(seq_pos);
                }

                // Go to next operation
                ml_pos++;
            }

            // Treat remaining unmodified cytosines as unmethylated
            while (seq_pos != END_OF_SEQUENCE) {
                seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                // Do not call methylation as these cytosines are inconclusive
            }

            break;
        } else {
            std::vector<std::string> elements = split(modification, ",");
            ml_pos += elements.size();
        }
    }

}

#endif // METHYLATION_H_
