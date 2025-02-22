//  reference.h
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

#ifndef MODALIB_REFERENCE_H
#define MODALIB_REFERENCE_H

#include <cctype>

struct AlignmentRefSequence {
private:
    char* seq;
    int start_pos;
    int length;
public:
    AlignmentRefSequence(char* seq, int start_pos, int length) {
        this->seq = seq;
        this->start_pos = start_pos;
        this->length = length;
    }
    ~AlignmentRefSequence() = default;  // Don't delete `seq` pointer
    char get_base(int abs_pos) const {
        return *(this->seq + abs_pos - this->start_pos);
    }
    Context get_seq_context(
            char first_base,
            char previous_base,
            char base,
            char next_base,
            char last_base,
            int c_conversion,
            int abs_pos,
            int strand) const {
        int pos = abs_pos - this->start_pos;
        if (pos >= this->length) return UNKNOWN_CONTEXT;
        char ref_base = static_cast<char>(std::toupper(*(this->seq + pos)));

        if ((strand & 1) == 1) {  // Forward strand
            if ((ref_base == 'C') && ((base == 'C') || (c_conversion && (base == 'T')))) {
                if (pos + 1 == this->length) return UNKNOWN_CONTEXT;
                char next_ref_base = static_cast<char>(std::toupper(*(this->seq + pos + 1)));
                if ((next_ref_base == 'G') && (next_base == 'G')) return CpG_CONTEXT;
                else if ((next_ref_base != 'G') && (next_base != 'G')) {
                    if (pos + 2 == this->length) return UNKNOWN_CONTEXT;
                    char last_ref_base = static_cast<char>(std::toupper(*(this->seq + pos + 2)));
                    if ((last_ref_base == 'G') && (last_base == 'G')) return CHG_CONTEXT;
                    else if ((last_ref_base != 'G') && (last_base != 'G')) return CHH_CONTEXT;
                }
            }
        } else if ((strand & 1) == 0) {  // Reverse strand
            if ((ref_base == 'G') && ((base == 'G') || (c_conversion && (base == 'A')))) {
                if (pos == 0) return UNKNOWN_CONTEXT;
                char previous_ref_base = static_cast<char>(std::toupper(*(this->seq + pos - 1)));
                if ((previous_ref_base == 'C') && (previous_base == 'C')) return CpG_CONTEXT;
                else if ((previous_ref_base != 'C') && (previous_base != 'C')) {
                    if (pos < 2) return UNKNOWN_CONTEXT;
                    char first_ref_base = static_cast<char>(std::toupper(*(this->seq + pos - 2)));
                    if ((first_ref_base == 'C') && (first_base == 'C')) return CHG_CONTEXT;
                    else if ((first_ref_base != 'C') && (first_base != 'C')) return CHH_CONTEXT;
                }
            }
        }
        return UNKNOWN_CONTEXT;
    }
};

#endif // MODALIB_REFERENCE_H
