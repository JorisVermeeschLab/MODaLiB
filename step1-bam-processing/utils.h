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

#include <string>

bool all_chars_the_same(const std::string& s) {
    if (s.empty()) return false;
    char c = s[0];
    for (int i = 1; i < s.length(); i++) {
        if (s[i] != c)
            return false;
    }
    return true;
}

#endif //MODALIB_UTILS_H
