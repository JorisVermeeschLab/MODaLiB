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
#include "reference.h"
#include "utils.h"

// Minimum length of a microsatellite motif (e.g., "AT")
const int MIN_MOTIF_LENGTH = 2;

// Maximum length of a microsatellite motif (e.g., "ATCCA")
const int MAX_MOTIF_LENGTH = 5;

// Minimum number of repeats to be considered a microsatellite
const int MIN_NUM_REPEATS = 5;

struct TandemRepeat {
    bool is_valid = false;
    uint16_t start = 0;
    uint16_t end = 0;
    std::string motif;
    uint16_t n_repeats = 0;
    inline bool spans(uint16_t pos) const { return (start <= pos) && (pos <= end); }
};


TandemRepeat find_microsatellite(std::vector<char> seq) {

    TandemRepeat tr;

    int n = seq.size();

    int counts[MAX_MOTIF_LENGTH - MIN_MOTIF_LENGTH + 1];
    for (int k = MIN_MOTIF_LENGTH; k <= MAX_MOTIF_LENGTH; k++) counts[k - MIN_MOTIF_LENGTH] = k;
    int best_i = 0;
    int best_n_repeats = 1;
    int best_motif_length = 1;

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < MAX_MOTIF_LENGTH - MIN_MOTIF_LENGTH + 1; k++) {
            int motif_length = k + MIN_MOTIF_LENGTH;
            if (i - motif_length >= 0) {
                if (seq.at(i) == seq.at(i - motif_length)) {
                    counts[k]++;
                } else {
                    counts[k] = motif_length;
                }
            }
            int n_repeats = counts[k] / motif_length;
            if (n_repeats * motif_length > best_n_repeats * best_motif_length) {
                best_n_repeats = n_repeats;
                best_motif_length = motif_length;
                best_i = i;
            }
        }
    }

    if (best_n_repeats >= MIN_NUM_REPEATS) {
        tr.is_valid = true;
        tr.n_repeats = best_n_repeats;
        tr.start = best_i + 1 - best_motif_length * best_n_repeats;
        tr.end = best_i;
        tr.motif = std::string(&seq[tr.end + 1 - best_motif_length], best_motif_length);

        // Handle homopolymeric sequences
        if (all_chars_the_same(tr.motif)) {
            tr.motif = tr.motif[0];
            tr.n_repeats = tr.end - tr.start + 1;
        }
    }
    return tr;
}

#endif //MODALIB_VNTR_H
