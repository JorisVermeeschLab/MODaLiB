//  end_motifs.h
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

#ifndef END_MOTIFS_H
#define END_MOTIFS_H

#include <algorithm>
#include <cstdint>
#include <vector>

#include "base.h"
#include "strand.h"


char reverse_complement(char base);


struct EndMotif {
private:
    char first_base = 'N';
    char second_base = 'N';
    char third_base = 'N';
    char fourth_base = 'N';
public:
    EndMotif() = default;
    EndMotif(char a, char b, char c, char d, int strand) {
        if ((strand & 1) == 1) {
            first_base = a;
            second_base = b;
            third_base = c;
            fourth_base = d;
        } else {
            first_base = reverse_complement(a);
            second_base = reverse_complement(b);
            third_base = reverse_complement(c);
            fourth_base = reverse_complement(d);
        }
    }
    bool is_valid() const {
        if (first_base == 'N') return false;
        if (second_base == 'N') return false;
        if (third_base == 'N') return false;
        if (fourth_base == 'N') return false;
        return true;
    }
    static uint8_t base_as_uint8(char c) {
        switch (c) {
            case 'A':
                return 0;
            case 'T':
                return 1;
            case 'C':
                return 2;
            case 'G':
                return 3;
            default:
                return 0;
        }
    }
    uint8_t as_uint8() const {
        uint8_t motif = 0;
        motif += base_as_uint8(first_base) << (unsigned int) 6;
        motif += base_as_uint8(second_base) << (unsigned int) 4;
        motif += base_as_uint8(third_base) << (unsigned int) 2;
        motif += base_as_uint8(fourth_base);
        return motif;
    }
    std::string as_str() const {
        return {first_base, second_base, third_base, fourth_base};
    }
};

#endif //END_MOTIFS_H
