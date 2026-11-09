//  site.h
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

#ifndef MODALIB_SITE_H
#define MODALIB_SITE_H

#include <cmath>
#include <functional>
#include <cstdint>
#include <utility>

#include "../fragment.h"
#include "overlap.h"


class Region {
private:
    int region_id;
    uint8_t chr_id;
    uint32_t start;
    uint32_t end;
    std::string motif;
    size_t n_repeats;
public:
    Region(
            int region_id,
            uint8_t chr_id,
            uint32_t start,
            uint32_t end,
            std::string motif,
            size_t n_repeats
            ) : region_id(region_id), chr_id(chr_id), start(start), end(end), motif(std::move(motif)), n_repeats(n_repeats) {}
    inline Overlap compute_overlap(const Fragment& f) const {
        return Overlap(region_id, chr_id, start, end, f.chr_id, f.start_pos, f.end_pos);
    }
    inline uint32_t get_start() const { return start; }
    inline uint32_t get_end() const { return end; }
    inline std::string get_motif() const { return motif; }
    inline size_t get_n_repeats() const { return n_repeats; }
};

#endif // MODALIB_SITE_H
