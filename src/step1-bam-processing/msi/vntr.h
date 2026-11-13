//  vntr.h
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

#ifndef MODALIB_VNTR_H
#define MODALIB_VNTR_H

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "../base.h"
#include "../reference.h"


// Minimum length of a microsatellite motif (e.g., "A")
const int MIN_MOTIF_LENGTH = 1;

// Maximum length of a microsatellite motif (e.g., "ATCCA")
const int MAX_MOTIF_LENGTH = 6;

// Minimum number of complete motif copies to be considered a microsatellite
const int MIN_NUM_REPEATS = 5;

// Minimum number of bases to be considered a microsatellite
const int MS_MIN_NUM_BASES = 15;

// Permit terminal partial motifs when there are at least four complete copies
// and the total motif-consistent run is at least 14 bases (e.g., ACT*4 + AC).
// This broadens the reference catalog: validate specificity for your use case.
const int MIN_FULL_COPIES_WITH_PARTIAL = 4;
const int MIN_PARTIAL_REPEAT_BASES = 14;


struct TandemRepeat {
    bool is_valid = false;
    // Reference coordinates for read-based calls: 1-based inclusive.
    // Sequence-local coordinates for find_microsatellite_in_reference(): 0-based inclusive.
    uint32_t start = 0;
    uint32_t end = 0;
    std::string motif;
    uint16_t n_repeats = 0;  // Complete motif copies
    uint8_t terminal_partial_bases = 0;
    int16_t indels = 0;
    inline bool spans(uint32_t pos) const {
        return (start <= pos + 1) && (pos <= end + 1);
    }
};


TandemRepeat find_microsatellite_in_reference(const char* seq, size_t n);


void find_microsatellites_in_reference(
        const char* seq,
        size_t seq_end,
        size_t seq_offset,
        std::vector<TandemRepeat>& mss);


void find_microsatellites(
        const std::vector<std::pair<uint32_t, BaseInfo>>& bases,
        std::vector<TandemRepeat>& mss);


void find_microsatellites(
        const std::vector<std::pair<uint32_t, BaseInfo>>& b1,
        const std::vector<std::pair<uint32_t, BaseInfo>>& b2,
        std::vector<TandemRepeat>& mss
        );

#endif //MODALIB_VNTR_H
