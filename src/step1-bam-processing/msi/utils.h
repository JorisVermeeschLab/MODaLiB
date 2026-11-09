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

#include <cctype>
#include <string>
#include <utility>


std::pair<std::string, size_t> parse_mss_motif(const std::string& s) {
    std::string motif;
    std::string n_repeats;
    for (char c: s) {
        if (isalpha(c)) {
            motif += c;
        } else if (isdigit(c)) {
            n_repeats += c;
        }
    }
    size_t i = 0;
    if (!n_repeats.empty()) {
        i = std::stoi(n_repeats);
    }
    return {motif, i};
}

#endif //MODALIB_UTILS_H
