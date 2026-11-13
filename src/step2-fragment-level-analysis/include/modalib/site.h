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
#include <cstdint>

#include "fragment_info.h"
#include "overlap.h"

namespace modalib {

class Region {
private:
    struct Coordinate {
        uint8_t chr;
        uint32_t pos;
    };

    int region_id;
    Coordinate start;
    Coordinate end;
public:
    Region(int region_id, uint8_t chr_id, uint32_t start, uint32_t end)
        : region_id(region_id), start{chr_id, start}, end{chr_id, end} {}
    int cmp(const FragmentInfo &r);
    inline Overlap compute_overlap(const FragmentInfo& f) const {
        return Overlap(region_id, start.chr, start.pos, end.pos, f.chr_id, f.start_pos, f.end_pos);
    }
    inline uint8_t get_chr_id() const { return start.chr; }
    inline uint32_t get_start() const { return start.pos; }
    inline uint32_t get_end() const { return end.pos; }
};

} // namespace modalib

#endif // MODALIB_SITE_H
