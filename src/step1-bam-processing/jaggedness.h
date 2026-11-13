//  jaggedness.h
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

#ifndef JAGGEDNESS_H_
#define JAGGEDNESS_H_

#include <fstream>
#include <string>

static constexpr size_t MAX_FRAGMENT_LENGTH = 500;

long overall_methylated_counts[MAX_FRAGMENT_LENGTH + 1][MAX_FRAGMENT_LENGTH + 1]{};
long overall_unmethylated_counts[MAX_FRAGMENT_LENGTH + 1][MAX_FRAGMENT_LENGTH + 1]{};

void increment_methylation_count(int base_position, int fragment_length, bool methylation_status) {
    if (base_position < 0) return;
    if ((base_position >= fragment_length) || (fragment_length > MAX_FRAGMENT_LENGTH)) return;
    if (methylation_status) {
        overall_methylated_counts[fragment_length][base_position]++;
    } else {
        overall_unmethylated_counts[fragment_length][base_position]++;
    }
}

void write_methylation_counts_to_file(const std::string& path, char sep = ',') {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Failed to open " + path);
    out << "fragment_length" << sep << "base_position" << sep << "methylated" << sep << "unmethylated" << '\n';
    for (std::size_t i = 0; i < MAX_FRAGMENT_LENGTH + 1; ++i) {
        for (std::size_t j = 0; j < MAX_FRAGMENT_LENGTH + 1; ++j) {
            if ((overall_methylated_counts[i][j] > 0) || (overall_unmethylated_counts[i][j] > 0)) {
                out << i << sep << j << sep << overall_methylated_counts[i][j] << sep << overall_unmethylated_counts[i][j] << '\n';
            }
        }
    }
}

#endif // JAGGEDNESS_H_
