//  utils.h
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

#ifndef MODALIB_UTILS_H
#define MODALIB_UTILS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

namespace modalib {

// Constant indicating that EOF has been reached
    const long END_OF_SAMPLE = std::numeric_limits<long>::max();


/**
 * Combine genomic coordinates to a single integer.
 *
 * @param chr_id Chromosome identifier
 * @param pos Position on chromosome (in bp).
 * @return Combined genomic coordinates.
 */
    inline long combine_coords(uint8_t chr_id, long pos) {
        if (pos == END_OF_SAMPLE) return END_OF_SAMPLE;
        else return static_cast<long>(chr_id) * 1e+9 + pos;
    }

/**
 * Parse chromosome name and return chromosome identifier.
 *
 * @param c Chromosome name.
 * @return Chromosome identifier.
 */
    inline uint8_t parse_chr_name(const char* c) {
        size_t n = strlen(c);
        uint8_t chr;

        // If string starts with "chr", ignore the 3 first characters
        if ((n > 3) && (c[0] == 'c') && (c[1] == 'h') && (c[2] == 'r')) {
            c += 3;
            n -= 3;
        }

        char *p;
        long converted = strtol(c, &p, 10);
        if (*p) { // If not numeric
            if (n == 1) {
                if (c[0] == 'X') chr = 23;
                else if (c[0] == 'Y') chr = 24;
                else if (c[0] == 'M') chr = 25;
                else chr = 0;
            } else if ((n == 3) && (strncmp(c, "EBV", 3) == 0)) {
                chr = 26;
            } else {
                chr = 0;
            }
        }
        else { // If numeric
            chr = (uint8_t) converted;
        }
        return chr;
    }

} // namespace modalib

#endif // MODALIB_UTILS_H
