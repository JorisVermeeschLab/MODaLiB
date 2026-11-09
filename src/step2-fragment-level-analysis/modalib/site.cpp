//  site.cpp
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

#include "modalib/site.h"
#include "modalib/fragment_info.h"


namespace modalib {

Region::Region(uint8_t chr, uint32_t start, uint32_t end) {
    this->start.chr = chr;
    this->start.pos = start;
    this->end.chr = chr;
    this->end.pos = end;
}

/**
 * Compare the position of a fragment with respect to the closest region.
 *
 * If the read has no intersection with the region, and is located left
 * to the region, 1 will be returned.
 * (r.end < start)
 * 
 * If the read has no intersection with the region, and is located right
 * to the region, -1 will be returned.
 * (r.start > end)
 * 
 * If the read starts in the region, -2 will be returned.
 * (r.start >= start) && (r.start <= end) && (r.end > end)
 * 
 * If the read ends in the region, 2 will be returned.
 * (r.start < start) && (r.end >= start) && (r.end <= end)
 *
 * If the read both starts and ends in the region, 3 will be returned.
 * (r.start >= start) && (r.end <= end)
 *
 * If the read spans the entire region, 4 will be returned.
 * (r.start < start) && (r.end > end)
 * 
 * @param r Read summary.
 * @return Status indicating the relative position of the read.
 */
int Region::cmp(const FragmentInfo &r) {
    uint8_t chr_id = this->start.chr;
    uint32_t start = this->start.pos;
    uint32_t end = this->end.pos;
    uint32_t r_end = r.end_pos;

    if (r.start_chr_id == chr_id) {
        if (r_end < start) {
            return 1;
        } else if (r.start_pos > end) {
            return -1;
        } else if ((r.start_pos >= start) && (r.start_pos <= end) && (r_end > end)) {
            return -2;
        } else if ((r.start_pos < start) && (r_end >= start) && (r_end <= end)) {
            return 2;
        } else if ((r.start_pos >= start) && (r_end <= end)) {
            return 3;
        } else if ((r.start_pos < start) && (r_end > end)) {
            return 4;
        } else {
            return 5;
        }
    } else if (r.start_chr_id < chr_id) {
        return 1;
    } else {
        return -1;
    }
}

} // namespace modalib