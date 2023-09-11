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

#include "reference.h"

// Minimum length of a microsatellite motif (e.g., "AT")
const int MIN_MOTIF_LENGTH = 2;

// Maximum length of a microsatellite motif (e.g., "ATCCA")
const int MAX_MOTIF_LENGTH = 5;

struct TandemRepeat {

};

TandemRepeat find_microsatellite(AlignmentRefSequence& seq) {

    int n = seq.get_length();

    int counts[MAX_MOTIF_LENGTH - MIN_MOTIF_LENGTH + 1];
    for (int k = MIN_MOTIF_LENGTH; k <= MAX_MOTIF_LENGTH; k++) counts[k - MIN_MOTIF_LENGTH] = k;
    int best_i = 0;
    int best_n_repeats = 1;
    int best_motif_length = 1;

    int start_abs_pos = seq.get_start_pos();

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < MAX_MOTIF_LENGTH - MIN_MOTIF_LENGTH + 1; k++) {
            int motif_length = k + MIN_MOTIF_LENGTH;
            if (i - motif_length >= 0) {
                if (seq.get_base(i + start_abs_pos) == seq.get_base(i + start_abs_pos - motif_length)) {
                    counts[k]++;
                } else {
                    counts[k] = motif_length;
                }
            }
            int n_repeats = counts[k] / motif_length;
            if (n_repeats >= 5) {
                if (n_repeats * motif_length > best_n_repeats * best_motif_length) {
                    best_n_repeats = n_repeats;
                    best_motif_length = motif_length;
                    best_i = i;
                }
            }
        }
    }
}

#endif //MODALIB_VNTR_H
