//  utils.cpp
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

#include "utils.h"

#include <cstdint>
#include <vector>


std::pair<uint8_t, uint8_t> methylation_counts_as_uint8_t(uint16_t n_methylated, uint16_t n_unmethylated) {
    if ((n_methylated == 0) && (n_unmethylated == 0)) return {0, 0};
    float meth_ratio = static_cast<float>(n_methylated) / static_cast<float>(n_methylated + n_unmethylated);
    if ((n_methylated > 255) && (n_methylated > n_unmethylated)) {
        return {255, static_cast<uint8_t>((1.0f - meth_ratio) * 255.0f / meth_ratio)};
    } else if ((n_unmethylated > 255) && (n_unmethylated > n_methylated)) {
        return {static_cast<uint8_t>(meth_ratio * 255.0f / (1.0f - meth_ratio)), 255};
    } else {
        return {n_methylated, n_unmethylated};
    }
}
