//  features.h
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

#ifndef _FEATURES_H__
#define _FEATURES_H__

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <valarray>
#include <vector>
#include <cstdlib>

#include "array_utils.h"
#include "fragment_info.h"
#include "bounds.h"
#include "sample.h"
#include "utils.h"

namespace modalib {

class FeatureMatrix {
private:
    size_t n;
    size_t m;
    bool normalize;
    std::valarray<long double> data;
    std::valarray<long double> counts;
public:
	FeatureMatrix() : n(0), m(0), normalize(false), data(0), counts(0) {}
	FeatureMatrix(size_t n, size_t m, bool normalize): n(n), m(m), normalize(normalize), data(0.0, n * m), counts(0.0, n * m) {}
	FeatureMatrix(size_t n, size_t m): FeatureMatrix(n, m, true) {}
	FeatureMatrix(const FeatureMatrix& other): n(other.n), m(other.m), normalize(other.normalize), data(other.data), counts(other.counts) {}
	double operator()(int i, int j) {
		long double count = this->counts[i * this->m + j];
		if (normalize) {
            if (count == 0) return 0.0;
            return this->data[i * this->m + j] / count;
        } else {
            return count;
        }
	}
	inline void add(float value, int i, int j, float weight) {
	    this->counts[i * this->m + j] += weight;
	    this->data[i * this->m + j] += value * weight;
	}
	inline void add_unweighted(float value, int i, int j, float weight) {
	    this->counts[i * this->m + j] += weight;
	    this->data[i * this->m + j] += value;
	}
	size_t n_rows() { return this->n; }
	size_t n_cols() { return this->m; }
};

std::unordered_map<uint8_t, std::vector<FeatureMatrix>> extract_all_features(
        std::vector<_Sample>& samples,
        std::unordered_map<uint8_t, std::vector<long>>& nuc_pos
        ) {

	size_t n = samples.size();

	std::unordered_map<uint8_t, std::vector<FeatureMatrix>> matrices;
	for (uint8_t chr_id = 0; chr_id < CHROMOSOME_BOUNDS.size(); ++chr_id) {
        size_t m = CHROMOSOME_BOUNDS.at(chr_id).second;
        matrices[chr_id] = std::vector<FeatureMatrix>({
            FeatureMatrix(n, m),  // Nucleosomes
            FeatureMatrix(n, m),  // Long-short ratios
            FeatureMatrix(n, m),  // Methylation
            FeatureMatrix(n, m, false),  // Frag counts
            FeatureMatrix(n, m)  // Mutation burden
        });
	}

    // TODO: reset samples

	Oscillator osc(nuc_pos);

    for (size_t i = 0; i < n; ++i) {

        _Sample& sample = samples.at(i);

        long current_combined_coord = 0;

    	while (!sample.eof()) {    		

            const FragmentInfo& frag_i = sample.get();

    		const uint8_t chr_id = frag_i.start_chr_id;

            // Get statistics about current frag
            const long coord_i = frag_i.start_pos;
            const long combined_coord_i = combine_coords(chr_id, coord_i);
            const long l_i = frag_i.size;
            const float depth_i = static_cast<float>(frag_i.n_sites());
            const float m_i = frag_i.get_meth_ratio();

            // Ensure fragments are sorted
            if (combined_coord_i < current_combined_coord) {
                std::cout << "Fragments are not properly sorted: (" << combined_coord_i << " < ";
                std::cout << current_combined_coord << ")" << std::endl;
            }

            // Get bin index
		    const size_t j = coord_i / 50000;

            // Update nucleosome features
		    if (j < matrices.at(chr_id).at(0).n_cols()) {
                float score = osc.phase_shift_score_wrt_reference(chr_id, static_cast<unsigned long>(coord_i));
                score += osc.phase_shift_score_wrt_reference(chr_id, static_cast<unsigned long>(coord_i + l_i));
                matrices.at(chr_id).at(0).add(0.5f * score, i, j, 1.0f);
            }

            // Update long short ratios
            if (j < matrices.at(chr_id).at(1).n_cols()) {
                if (l_i > 150) {
                    matrices.at(chr_id).at(1).add_unweighted(1.0f, i, j, 0.0f);
                } else {
                    matrices.at(chr_id).at(1).add_unweighted(0.0f, i, j, 1.0f);
                }
            }

            // Update genome-wide methylation features
            if (j < matrices.at(chr_id).at(2).n_cols()) {
                matrices.at(chr_id).at(2).add(m_i, i, j, depth_i);
            }

            // Update frag counts
            if (j < matrices.at(chr_id).at(3).n_cols()) {
                matrices.at(chr_id).at(3).add(1.0f, i, j, 1.0f);
            }

            // Update mutation burden features
            if (j < matrices.at(chr_id).at(4).n_cols()) {
                matrices.at(chr_id).at(4).add(frag_i.get_mutation_burden(), i, j, frag_i.s_count);
            }

    		// Move to next fragment
    		sample.next();

    	}
    }


	return matrices;
}

} // namespace modalib

#endif // _FEATURES_H__
