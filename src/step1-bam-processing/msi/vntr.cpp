//  vntr.cpp
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

#include "vntr.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "../base.h"
#include "../reference.h"


TandemRepeat find_microsatellite_in_reference(const char* seq, size_t n) {

    TandemRepeat tr;

    if (n == 0) return tr;

    int counts[MAX_MOTIF_LENGTH - MIN_MOTIF_LENGTH + 1];
    for (int k = MIN_MOTIF_LENGTH; k <= MAX_MOTIF_LENGTH; k++) counts[k - MIN_MOTIF_LENGTH] = k;

    for (int i = 0; i < n; i++) {
        // for (int k = 0; k < MAX_MOTIF_LENGTH - MIN_MOTIF_LENGTH + 1; k++) {
        for (int k = MAX_MOTIF_LENGTH - MIN_MOTIF_LENGTH; k >= 0; k--) {
            int motif_length = k + MIN_MOTIF_LENGTH;
            if (i - motif_length >= 0) {
                if ((seq[i] != seq[i - motif_length]) || (seq[i] == 'N')) {
                    counts[k] = motif_length;

                    if (tr.is_valid && (tr.motif.length() == motif_length)) {
                        return tr;
                    }
                } else {
                    counts[k]++;
                }
            }
            if (counts[k] >= MS_MIN_NUM_BASES) {
                int n_repeats = counts[k] / motif_length;
                tr.is_valid = true;
                tr.n_repeats = n_repeats;
                tr.start = std::max(i + 1 - counts[k], 0);
                tr.end = i;
                tr.motif = std::string(&seq[tr.start], motif_length);

                // Ensure that microsatellite is entirely contained in the sequence
                if (tr.end >= n - 1) tr.is_valid = false;
                else if (seq[tr.end + 1] == 'N') tr.is_valid = false;
                if (tr.start <= 0) tr.is_valid = false;
                else if (seq[tr.start - 1] == 'N') tr.is_valid = false;
            }
        }
    }
    return tr;
}


void find_microsatellites_in_reference(const char* seq, size_t seq_end, size_t seq_offset, std::vector<TandemRepeat>& mss) {
    size_t seq_start = 0;
    TandemRepeat tr;
    do {
        tr = find_microsatellite_in_reference(seq + seq_start, seq_end - seq_start);
        if (tr.is_valid) {
            size_t new_seq_start = seq_start + tr.end + 1;
            tr.start += seq_start + seq_offset;
            tr.end += seq_start + seq_offset;
            mss.push_back(tr);
            seq_start = new_seq_start;
            if (seq_start >= seq_end) break;
        }
    } while (tr.is_valid);
}


void find_microsatellites(const std::vector<std::pair<uint32_t, BaseInfo>>& bases, std::vector<TandemRepeat>& mss) {
    if (bases.empty()) return;
    char* seq = new char[bases.size()];
    size_t seq_offset = bases.front().first;
    size_t seq_end = 0;
    for (const auto& pair : bases) {
        if (pair.second.deleted) continue;
        if (pair.second.soft_clipped) seq[seq_end++] = 'N';
        else {
            seq[seq_end++] = pair.second.c;
        }
    }
    find_microsatellites_in_reference(seq, seq_end, seq_offset, mss);
    delete[] seq;
}


void find_microsatellites(
        const std::vector<std::pair<uint32_t, BaseInfo>>& b1,
        const std::vector<std::pair<uint32_t, BaseInfo>>& b2,
        std::vector<TandemRepeat>& mss
        ) {
    if (b1.size() + b2.size() == 0) return;
    char* seq = new char[b1.size() + b2.size()];
    size_t seq_offset;
    if (!b1.empty()) seq_offset = b1.front().first;
    else seq_offset = b2.front().first;
    size_t seq_end = 0;
    uint32_t last = 0;
    for (const auto& pair : b1) {
        last = pair.first;
        if (pair.second.deleted) continue;
        if (pair.second.soft_clipped) seq[seq_end++] = 'N';
        else seq[seq_end++] = pair.second.c;
    }
    for (const auto& pair : b2) {
        if (last < pair.first) {
            last = pair.first;
            if (pair.second.deleted) continue;
            if (pair.second.soft_clipped) seq[seq_end++] = 'N';
            else seq[seq_end++] = pair.second.c;
        }
    }
    find_microsatellites_in_reference(seq, seq_end, seq_offset, mss);
    delete[] seq;
}
