//  nucleosomes.h
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

#ifndef MODALIB_NUCLEOSOMES_H
#define MODALIB_NUCLEOSOMES_H

#include "bounds.h"
#include "utils.h"
#include "fragment_info.h"

#include <valarray>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace modalib {

    // https://download.cncb.ac.cn/nucmap/organisms/v1/Homo_sapiens/byDataType/Nucleosome_peaks_DANPOS/
    // https://academic.oup.com/nar/article/47/D1/D163/5134334

class RefNucMap {
private:

    // Nucleosome positions per chromosome.
    // Keys are chromosome ids, while values are vectors of called nucleosome positions.
    // The size of each vector corresponds to the number of nucleosomes called on each chromosome.
    std::unordered_map<uint8_t, std::vector<long>> nuc_pos;

    // Angles w.r.t. nucleosome positions for each chromosome.
    // Size of each vector is the size of the corresponding chromosome (in bp).
    std::unordered_map<uint8_t, std::vector<float>> chr_angles;

public:
    RefNucMap() = default;
    RefNucMap(const std::string& filepath) {

        long last_combined_coords = 0;

        // Parse input bed file
        std::ifstream file(filepath);
        std::vector<std::string> line_elements;
        for (std::string line; std::getline(file, line);) {
            line_elements.clear();
            std::istringstream ss(line);
            for (std::string value; std::getline(ss, value, '\t');) {
                line_elements.push_back(std::move(value));
            }
            if (line_elements.size() == 4) {

                // Parse nucleosome identifier
                uint8_t chr_id = parse_chr_name(line_elements.at(0).c_str());

                // Skip nucleosome if unknown chromosome
                if (chr_id == 0) continue;

                // Get nucleosome center position
                uint32_t start = std::stoi(line_elements.at(1));
                uint32_t end = std::stoi(line_elements.at(2));
                uint32_t center = (start + end) / 2;

                // Make sure bed file is properly sorted
                long combined_coords = combine_coords(chr_id, center);
                if (combined_coords < last_combined_coords) {
                    throw std::domain_error(
                            "Error while parsing reference nucleosome map. "
                            "Input bed file should be sorted (chr1, chr2, ..., chr11, ..., chr22, chrX, chrY, chrM)");
                }
                last_combined_coords = combined_coords;

                // Store position center position
                if (nuc_pos.find(chr_id) == nuc_pos.end()) {
                    nuc_pos[chr_id] = std::vector<long>();
                }
                nuc_pos[chr_id].push_back(center);

                // Allocate memory for nucleosome scores
                if (chr_angles.find(chr_id) == chr_angles.end()) {
                    if (chr_id < CHROMOSOME_BOUNDS.size()) {
                        long chr_size = CHROMOSOME_BOUNDS.at(chr_id).second * BIN_SIZE;
                        chr_angles[chr_id] = std::vector<float>(chr_size, 0.0f);
                    }
                }
            }
        }

        // Update nucleosome scores
        for (auto& it: chr_angles) {
            uint8_t chr_id = it.first;
            std::vector<float> &ca = it.second;
            std::vector<long> &pos = nuc_pos.at(chr_id);
            size_t n = pos.size();
            if (n <= 1) continue;
            for (int i = 1; i < n; i++) {
                const long x1 = pos[i - 1];
                const long x2 = pos[i];
                if (x2 >= ca.size()) continue;
                if ((x2 - x1 > 500) || (x2 - x1 < 100)) {
                    for (long j = 0; j < x2 - x1; j++) {
                        ca[x1 + j] = 0.0f;
                    }
                }
                for (long j = 0; j < x2 - x1; j++) {
                    ca[x1 + j] = 2.0f * static_cast<float>(j) / static_cast<float>(x2 - x1);
                }
            }
        }
        std::cout << "Reference nucleosome map loaded." << std::endl;
    }

    inline float angle(uint8_t chr_id, unsigned long genomic_coord) const {
        const std::vector<float>& ca = this->chr_angles.at(chr_id);
        return (genomic_coord >= ca.size()) ? 0.0f : ca[genomic_coord];
    }

    std::valarray<float> phase_shift_scores(
            uint8_t chr_id, const std::valarray<FragmentInfo>& frags, const std::vector<int>& idx, const FragmentInfo& frag_i
    ) const {

        std::valarray<float> pss(0.0f, idx.size());

        if (this->chr_angles.find(chr_id) == this->chr_angles.end()) return pss;
        const std::vector<float>& ca = this->chr_angles.at(chr_id);
        if (frag_i.end_pos >= ca.size()) return pss;

        const float phi2 = ca.at(frag_i.start_pos);
        const float psi2 = ca.at(frag_i.end_pos);

        for (size_t k = 0; k < idx.size(); k++) {

            const int j = idx[k];

            if (frags[j].end_pos >= ca.size()) continue; // TODO: set weight to 0 (add weights as function arg)

            // Compute score for start positions
            float phi = std::fabs(ca.at(frags[j].start_pos) - phi2);
            phi = 1.0f - std::fmin(phi, 2.0f - phi);

            // Compute score for end positions
            float psi = std::fabs(ca.at(frags[j].end_pos) - psi2);
            psi = 1.0f - std::fmin(psi, 2.0f - psi);

            // Compute total score
            pss[k] = 0.5 * phi + 0.5 * psi;
        }
        return pss;
    }

    float phase_shift_score_wrt_reference(uint8_t chr_id, unsigned long x) const {
        float phi = this->angle(chr_id, x);
        phi = std::fmin(phi, 2.0f - phi);
        return 1.0f - phi;
    }
};

} // namespace modalib

#endif // MODALIB_NUCLEOSOMES_H
