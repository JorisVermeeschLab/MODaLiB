//  base.h
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

#ifndef MODALIB_BASE_H
#define MODALIB_BASE_H

#include <cstdint>

const uint8_t BASE_TYPE_NORMAL = 0;
const uint8_t BASE_TYPE_PUTATIVE_SNP = 1;
const uint8_t BASE_TYPE_METHYLATED_CPG = 2;
const uint8_t BASE_TYPE_UNMETHYLATED_CPG = 3;
const uint8_t BASE_TYPE_UNINFORMATIVE = 4;
const uint8_t BASE_TYPE_INSERTED = 5;
const uint8_t BASE_TYPE_DELETED = 6;
const uint8_t BASE_TYPE_SOFT_CLIPPED = 7;
const uint8_t BASE_TYPE_REF_SKIPPED = 8;
uint8_t BASE_TYPE_AS_CHARS[9] = {' ', 'X', 'M', 'U', '0', '/', '\\', 'S', 'R'};

struct BaseInfo {

    // Char representation
    char c = 'N';

    // Phred score
    uint8_t phred_score = 0;

    // Base type (normal, mutated, etc.)
    uint8_t type = BASE_TYPE_NORMAL;

    inline bool is_c_or_g() const { return (c == 'C') or (c == 'G'); }
    inline char type_as_char() const { return BASE_TYPE_AS_CHARS[type]; }
};

#endif //MODALIB_BASE_H
