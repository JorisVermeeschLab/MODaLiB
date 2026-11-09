//  overlap.h
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

#ifndef MODALIB_OVERLAP_H
#define MODALIB_OVERLAP_H

#include <cstdint>

namespace modalib {


const int RPOS_LEFT_NO_OVERLAP = 1;
const int RPOS_RIGHT_NO_OVERLAP = -1;
const int RPOS_LEFT_END_OUTSIDE = 2;
const int RPOS_RIGHT_END_OUTSIDE = -2;
const int RPOS_INSIDE = 3;
const int RPOS_FLANKING_ENDS = 4;
const int RPOS_UNDEFINED = 5;


class Overlap {
public:

    // Index of the region being overlapped.
    int region_id;

    // How the fragment is positioned with respect to the region.
    // e.g., with the left end outside the region and the right end inside.
    int rpos_type;

    // Whether the left end of the fragment is located left (outside) to the region,
    // indicating that the fragment potentially overlaps the previous region.
    bool may_overlap_previous_region;

    // Whether the right end of the fragment is located right (outside) to the region,
    // indicating that the fragment potentially overlaps the next region.
    bool may_overlap_next_region;

    // Overlap fraction, expressed in terms of the fragment.
    // For example, if both ends of a fragment are located inside the region,
    // then overlap is 1.
    float fraction;

    // Results of the region search
    bool found_contig = true;
    bool reached_contig_end = false;

    Overlap() : region_id(0), rpos_type(RPOS_UNDEFINED), may_overlap_previous_region(false), may_overlap_next_region(false), fraction(0.0f) {}
    Overlap(int region_id, uint8_t chr_id, uint32_t start, uint32_t end, uint8_t r_chr_id, uint32_t r_start, uint32_t r_end);
    inline bool left_end_inside() const { return (rpos_type == RPOS_RIGHT_END_OUTSIDE) || (rpos_type == RPOS_INSIDE); }
    inline bool right_end_inside() const { return (rpos_type == RPOS_LEFT_END_OUTSIDE) || (rpos_type == RPOS_INSIDE); }
};

} // namespace modalib

#endif // MODALIB_OVERLAP_H
