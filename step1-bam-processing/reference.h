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

#ifndef REFERENCE_H_
#define REFERENCE_H_

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
    char get_base(int ref_pos) {
        return *(this->seq + ref_pos - this->start_pos);
    }
    int is_cpg(int ref_pos) {
        int pos = ref_pos - this->start_pos;
        if (pos >= this->length) return 0;
        char base = static_cast<char>(std::toupper(*(this->seq + pos)));
        if (base == 'C') {
            if (pos + 1 == this->length) return 0;
            char next_base = static_cast<char>(std::toupper(*(this->seq + pos + 1)));
            if (next_base == 'G') return 1;
            return 0;
        } else if (base == 'G') {
            if (pos == 0) return 0;
            char previous_base = static_cast<char>(std::toupper(*(this->seq + pos - 1)));
            if (previous_base == 'C') return -1;
            return 0;
        }
        return 0;
    }
};

#endif // READ_INFO_H_
