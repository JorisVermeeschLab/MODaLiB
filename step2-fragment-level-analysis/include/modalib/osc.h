//  osc.h
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

#ifndef _OSC_H__
#define _OSC_H__

#include "modalib/fragment_info.h"

#include <fstream>
#include <iostream>
#include <string>
#include <cstdint>
#include <utility>
#include <valarray>
#include <vector>
#include <unordered_map>

namespace modalib {

class Oscillator {

private:

    // Nucleosome positions per chromosome.
    // Keys are chromosome ids, while values are vectors of called nucleosome positions.
    // The size of each vector corresponds to the number of nucleosomes called on each chromosome.
    std::unordered_map<uint8_t, std::vector<long>> nuc_pos;

    // Angles w.r.t. nucleosome positions for the each chromosome.
    // Size of each vector is the size of the corresponding chromosome (in bp).
    std::unordered_map<uint8_t, std::vector<float>> chr_angles;

public:

    // No default constructor (to make sure it is properly initialized)
    Oscillator() = delete;


	Oscillator(std::unordered_map<uint8_t, std::vector<long>>& nuc_pos): nuc_pos(), chr_angles(0) {
	    for (auto& it : nuc_pos) {

            const uint8_t chr_id = it.first;

	        std::cout << "Loading reference nucleosome map for chr" << 1 * chr_id << "... ";

            const std::vector<long>& pos = it.second;
	        this->nuc_pos[chr_id] = pos;

            const int n = pos.size();
            const int size = (n == 0) ? 0 : pos.at(n - 1);
            this->chr_angles[chr_id] = std::vector<float>(size);
            std::vector<float>& ca = this->chr_angles[chr_id];
            if (n == 0) continue;
            for (int i = 1; i < n; i++) {
                const long x1 = pos[i - 1];
                const long x2 = pos[i];
                if ((x2 - x1 > 500) || (x2 - x1 < 100)) {
                    for (long j = 0; j < x2 - x1; j++) {
                        ca[x1 + j] = 0.0f;
                    }
                }
                for (long j = 0; j < x2 - x1; j++) {
                    ca[x1 + j] = 2.0f * static_cast<float>(j) / static_cast<float>(x2 - x1);
                }
            }

            std::cout << "done." << std::endl;
	    }
	}
	std::valarray<float> phase_shift_scores(
	    uint8_t chr_id, const std::valarray<FragmentInfo>& fragments, const std::vector<int>& idx, const FragmentInfo& frag_i
	) const;
	float phase_shift_score_wrt_reference(uint8_t chr_id, unsigned long x) const;
    inline float angle(uint8_t chr_id, unsigned long genomic_coord) const {
        const std::vector<float>& ca = this->chr_angles.at(chr_id);
        return (genomic_coord >= ca.size()) ? 0.0 : ca[genomic_coord];
    }
};

} // namespace modalib

#endif // _OSC_H__
