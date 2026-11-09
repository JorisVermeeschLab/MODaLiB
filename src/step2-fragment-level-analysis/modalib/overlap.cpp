//  overlap.cpp
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

#include <iostream>
#include <cmath>

#include "modalib/overlap.h"

namespace modalib {

/**
 * Compare the position of a fragment with respect to the closest region.
 *
 * If the read has no intersection with the region, and is located left
 * to the region, RPOS_LEFT_NO_OVERLAP will be returned.
 * (r.end < start)
 * 
 * If the read has no intersection with the region, and is located right
 * to the region, RPOS_RIGHT_NO_OVERLAP will be returned.
 * (r.start > end)
 * 
 * If the read starts in the region, RPOS_RIGHT_END_OUTSIDE will be returned.
 * (r.start >= start) && (r.start <= end) && (r.end > end)
 * 
 * If the read ends in the region, RPOS_LEFT_END_OUTSIDE will be returned.
 * (r.start < start) && (r.end >= start) && (r.end <= end)
 *
 * If the read both starts and ends in the region, RPOS_INSIDE will be returned.
 * (r.start >= start) && (r.end <= end)
 *
 * If the read spans the entire region, RPOS_FLANKING_ENDS will be returned.
 * (r.start < start) && (r.end > end)
 */
Overlap::Overlap(int region_id, uint8_t chr_id, uint32_t start, uint32_t end, uint8_t r_chr_id, uint32_t r_start, uint32_t r_end) {
    //std::cout << "CMP: " << 1*region_id << ", " << 1 * chr_id << ", " << start << ", " << end << ", " << 1* r_chr_id << ", " << r_start << ", " << r_end << std::endl;
    this->region_id = region_id;
    fraction = 0.0f;
    may_overlap_previous_region = true;
    may_overlap_next_region = true;
    if (r_chr_id == chr_id) {
        if (r_end < start) {
            // r_start r_end | |
            rpos_type = RPOS_LEFT_NO_OVERLAP;
            may_overlap_next_region = false;
        } else if (end < r_start) {
            // | | r_start r_end
            rpos_type = RPOS_RIGHT_NO_OVERLAP;
            may_overlap_previous_region = false;
        } else if ((start <= r_start) && (r_start <= end) && (end < r_end)) {
            // | r_start | r_end
            fraction = static_cast<float>(end + 1 - r_start) / static_cast<float>(r_end + 1 - r_start);
            rpos_type = RPOS_RIGHT_END_OUTSIDE;
            may_overlap_previous_region = false;
        } else if ((r_start < start) && (start <= r_end) && (r_end <= end)) {
            //  r_start | r_end |
            fraction = static_cast<float>(r_end + 1 - start) / static_cast<float>(r_end + 1 - r_start);
            rpos_type = RPOS_LEFT_END_OUTSIDE;
        } else if ((start <= r_start) && (r_end <= end)) {
            // | r_start r_end |
            fraction = 1.0f;
            rpos_type = RPOS_INSIDE;
            may_overlap_previous_region = false;
            may_overlap_next_region = false;
        } else if ((r_start < start) && (end < r_end)) {
            // r_start | | r_end
            fraction = static_cast<float>(end + 1 - start) / static_cast<float>(r_end + 1 - r_start);
            rpos_type = RPOS_FLANKING_ENDS;
            may_overlap_previous_region = true;
        } else {
            rpos_type = RPOS_UNDEFINED;
        }
    } else if (r_chr_id < chr_id) {
        rpos_type = RPOS_LEFT_NO_OVERLAP;
        may_overlap_next_region = false;
    } else {
        rpos_type = RPOS_RIGHT_NO_OVERLAP;
        may_overlap_previous_region = false;
    }

    if (std::isnan(fraction)) fraction = 0.0f;
}

} // namespace modalib
