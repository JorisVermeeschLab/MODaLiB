//  utils.h
//
//  Copyright 2024 Antoine Passemiers <antoine.passemiers@gmail.com>
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

#include <cstdint>
#include <utility>
#include <limits.h>


#if LONG_MAX > 4000000000
typedef long genome_size_t;
#else
typedef long long genome_size_t;
#endif


inline genome_size_t combine_chromosomal_coordinates(uint8_t chr_id, uint32_t pos) {
    return static_cast<genome_size_t>(chr_id) * 1000000000L + static_cast<genome_size_t>(pos);
}


inline std::pair<uint8_t, uint32_t> uncombine_chromosomal_coordinates(genome_size_t coords) {
    auto chr_id = static_cast<uint8_t>(coords / 1000000000L);
    auto pos = static_cast<uint32_t>(coords % 1000000000L);
    return {chr_id, pos};
}


#endif // MODALIB_UTILS_H
