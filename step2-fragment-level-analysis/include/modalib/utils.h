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

#ifndef UTILS_H_
#define UTILS_H_

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

} // namespace modalib

#endif // UTILS_H_
