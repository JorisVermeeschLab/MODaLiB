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

#include "htslib/sam.h"
#include "htslib/faidx.h"

#include <sstream>
#include <string>
#include <utility>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <utility>

#include "utils.h"


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

inline int get_base(uint8_t* seq, int pos, int seq_len, int strand) {
    if (!(strand & 1)) {
        pos = seq_len - 1 - pos;
    }
    if ((pos < 0) or (pos >= seq_len)) return -1;
    return seq_nt16_str[bam_seqi(seq, pos)];
}

inline bool is_cytosine(int base, int strand) {
    if (strand & 1) return base == 2;
    else return base == 4;
}

inline bool is_guanine(int base, int strand) {
    if (strand & 1) return base == 4;
    else return base == 2;
}

bool is_cpg(uint8_t* seq, int pos, int seq_len, int strand) {
    if (pos >= seq_len - 1) return false;
    int base = get_base(seq, pos, seq_len, strand);
    int next_base = get_base(seq, pos + 1, seq_len, strand);
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

std::pair<int, int> methylation_from_aux(bam1_t *b, uint8_t ml_cutoff) {

    int strand = get_strand(b);

    uint8_t* ml = bam_aux_get(b, "ML");
    char* mm = bam_aux2Z(bam_aux_get(b, "MM"));

    uint8_t* seq = bam_get_seq(b);
    int seq_len = b->core.l_qseq;

    std::pair<int, int> counts(0, 0);

    int seq_pos = -1;
    std::vector<std::string> modifications = split(mm, ";");

    int ml_pos = 0;

    for (std::string modification : modifications) {

        if (modification.rfind("C+h", 0) == 0) {

            std::vector<std::string> elements = split(modification, ",");

            for (int i = 1; i < elements.size(); i++) {

                for (int k = 0; k < std::stoi(elements.at(i)); k++) {
                    seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                    if (is_cpg(seq, seq_pos, seq_len, strand)) counts.second++;
                }

                seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                if (is_cpg(seq, seq_pos, seq_len, strand)) {
                    if (ml[ml_pos] >= ml_cutoff) counts.first++;
                    else counts.second++;
                }

                ml_pos++;
            }

            while (seq_pos != END_OF_SEQUENCE) {
                seq_pos = move_to_next_cytosine(seq, seq_pos, seq_len, strand);
                if (is_cpg(seq, seq_pos, seq_len, strand)) counts.second++;
            }

            break;
        } else {
            std::vector<std::string> elements = split(modification, ",");
            ml_pos += elements.size();
        }
    }

    return counts;
}


#endif // METHYLATION_H_
