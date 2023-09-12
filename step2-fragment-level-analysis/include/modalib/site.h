//  site.h
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

#ifndef MODALIB_SITE_H
#define MODALIB_SITE_H

#include <cmath>
#include <functional>
#include <stdint.h>

#include "fragment_info.h"

namespace modalib {

struct Site {
    uint8_t chr;
    uint32_t pos;
};

class Region {
public:
    Site start;
    Site end;
    Region(uint8_t chr, uint32_t start, uint32_t end);
    int cmp(const FragmentInfo&);
    bool operator==(const Region &other) const {
        bool v = (start.chr == other.start.chr);
        v &= (start.pos == other.start.pos);
        v &= (end.chr == other.end.chr);
        v &= (end.pos == other.end.pos);
        return v;
    }
};

struct RegionHashFunction {
    std::hash<uint8_t> hash_func8;
    std::hash<uint32_t> hash_func32;
    // TODO: optimize hash function: CRC32 hash algorithm?
    std::size_t operator()(const Region &r) const {
        std::size_t h = 157;  // TODO
        h = h ^ hash_func8(r.start.chr);
        h = h ^ hash_func32(r.start.pos);
        h = h ^ hash_func8(r.end.chr);
        h = h ^ hash_func32(r.end.pos);
        return h;
    }
};

} // namespace modalib

#endif // MODALIB_SITE_H
