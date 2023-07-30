//  osc.cpp
//
//  Copyright 2022 Antoine Passemiers <antoine.passemiers@gmail.com>
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

#include "modalib/osc.h"
#include "modalib/array_utils.h"


namespace modalib {

std::valarray<float> Oscillator::phase_shift_scores(
        uint8_t chr_id, const std::valarray<FragmentInfo>& frags, const std::vector<int>& idx, const FragmentInfo& frag_i
        ) const {

    std::valarray<float> pss(0.0f, idx.size());

    const std::vector<float>& ca = this->chr_angles.at(chr_id);
    if (frag_i.end_pos >= ca.size()) return pss;

    const float phi2 = ca.at(frag_i.start_pos);
    const float psi2 = ca.at(frag_i.end_pos);
    for (size_t k = 0; k < idx.size(); k++) {
        const int j = idx[k];

        if (frags[j].end_pos >= ca.size()) continue; // TODO: set weight to 0 (add weights as function arg)

        // Compute score for start positions
        float phi = std::fabs(ca[frags[j].start_pos] - phi2);
        phi = 1.0f - std::fmin(phi, 2.0f - phi);

        // Compute score for end positions
        float psi = std::fabs(ca[frags[j].end_pos] - psi2);
        psi = 1.0f - std::fmin(psi, 2.0f - psi);

        // Compute total score
        pss[k] = 0.5 * phi + 0.5 * psi;
    }
    return pss;
}

float Oscillator::phase_shift_score_wrt_reference(uint8_t chr_id, unsigned long x) const {
    float phi = this->angle(chr_id, x);
    phi = std::fmin(phi, 2.0f - phi);
    return 1.0f - phi;
}

} // namespace modalib
