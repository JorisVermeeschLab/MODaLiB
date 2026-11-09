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
#include "fragment_info.h"
#include "utils.h"

#include <valarray>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace modalib {

// https://download.cncb.ac.cn/nucmap/organisms/v1/Homo_sapiens/byDataType/Nucleosome_peaks_DANPOS/
// https://academic.oup.com/nar/article/47/D1/D163/5134334

class RefNucMap {
public:
    static constexpr float NULL_ANGLE = -1.0f;
    static constexpr float NULL_SCORE = -1.0f;

private:

    // Whether to return distance to closest nucleosome instead of phase score
    bool distance = false;

    // Nucleosome positions per chromosome.
    // Keys are chromosome ids, while values are vectors of called nucleosome positions.
    // The size of each vector corresponds to the number of nucleosomes called on each chromosome.
    std::unordered_map<uint8_t, std::vector<long>> nuc_pos;

    // Angles w.r.t. nucleosome positions for each chromosome.
    // Size of each vector is the size of the corresponding chromosome (in bp).
    std::unordered_map<uint8_t, std::vector<float>> chr_angles;

public:
    RefNucMap() = default;
    RefNucMap(const std::string& filepath): RefNucMap(filepath, false) {}
    RefNucMap(const std::string& filepath, bool distance) {
        this->distance = distance;

        genome_size_t last_combined_coords = 0;

        int n_nucleosomes = 0;

        // Parse input bed file
        std::ifstream file(filepath);
        std::vector<std::string> line_elements;
        for (std::string line; std::getline(file, line);) {

            line_elements.clear();
            std::istringstream ss(line);
            for (std::string value; std::getline(ss, value, '\t');) {
                line_elements.push_back(std::move(value));
            }
            if (line_elements.size() >= 2) {

                // Parse nucleosome identifier
                uint8_t chr_id = parse_chr_name(line_elements.at(0).c_str());

                // Skip nucleosome if unknown chromosome
                if (chr_id == 0) continue;

                // Get nucleosome center position
                long center = std::stoi(line_elements.at(1));

                // Make sure bed file is properly sorted
                genome_size_t combined_coords = combine_coords(chr_id, center);

                if (combined_coords < last_combined_coords) {
                    throw std::domain_error(
                            "Error while parsing reference nucleosome map. "
                            "Input bed file should be sorted (chr1, chr2, ..., chr11, ..., chr22, chrX, chrY)");
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
                        unsigned long chr_size = CHROMOSOME_BOUNDS.at(chr_id).second * BIN_SIZE;
                        chr_angles[chr_id] = std::vector<float>(chr_size, RefNucMap::NULL_ANGLE);
                    }
                }

                n_nucleosomes++;
            }
        }

        std::cout << "Found " << n_nucleosomes << " nucleosome(s) in reference map." << std::endl;

        // Update nucleosome scores
        for (auto& it: chr_angles) {
            uint8_t chr_id = it.first;
            std::vector<float> &ca = it.second;
            std::vector<long> &pos = nuc_pos.at(chr_id);
            size_t n = pos.size();
            if (n <= 1) continue;
            for (size_t i = 1; i < n; i++) {
                const unsigned long x1 = pos[i - 1];
                const unsigned long x2 = pos[i];
                if (x1 > x2) throw std::runtime_error("Nucleosome reference map file is not sorted!");
                if (x2 >= ca.size()) continue;
                if (x2 - x1 <= 500) {
                    // x1     ...   x2
                    // x1+0   ...   x1+x2-x1+1
                    //    0
                    for (unsigned long j = 0; j < x2 - x1; j++) {
                        if (!distance) {
                            ca[x1 + j] = 2.0f * static_cast<float>(j) / static_cast<float>(x2 - x1);
                        } else {
                            ca[x1 + j] = std::fmin(j, x2 - x1 - j);
                        }
                    }
                }
            }
        }
        std::cout << "Reference nucleosome map loaded." << std::endl;
    }

    inline float angle(uint8_t chr_id, unsigned long genomic_coord) const {
        if (this->chr_angles.find(chr_id) == this->chr_angles.end()) {
            return RefNucMap::NULL_ANGLE;
        } else {
            const std::vector<float>& ca = this->chr_angles.at(chr_id);
            return RefNucMap::angle(ca, genomic_coord);
        }
    }

    inline static float angle(const std::vector<float>& ca, unsigned long genomic_coord) {
        if (genomic_coord >= ca.size()) {
            return RefNucMap::NULL_ANGLE;
        } else {
            return ca[genomic_coord];
        }
    }

    inline static float default_angle(unsigned long genomic_coord) {
        return 2.0f * static_cast<float>(genomic_coord % 167) / 167.0f;
    }

    std::valarray<float> phase_shift_scores(
            uint8_t chr_id, const std::valarray<FragmentInfo>& frags, const std::vector<size_t>& idx, const FragmentInfo& frag_i
    ) const {

        std::valarray<float> pss(RefNucMap::NULL_SCORE, idx.size());

        // Ensure contig is available in reference map
        if (this->chr_angles.find(chr_id) == this->chr_angles.end()) return pss;
        const std::vector<float>& ca = this->chr_angles.at(chr_id);

        // Ensure positions are available in reference map
        const float phi2 = RefNucMap::angle(ca, frag_i.start_pos);
        const float psi2 = RefNucMap::angle(ca, frag_i.end_pos);

        for (size_t k = 0; k < idx.size(); k++) {

            const size_t j = idx[k];

            // Compute score for start positions
            float phi = std::fabs(RefNucMap::angle(ca, frags[j].start_pos) - phi2);
            if (phi == RefNucMap::NULL_ANGLE) continue;
            phi = 1.0f - std::fmin(phi, 2.0f - phi);

            // Compute score for end positions
            float psi = std::fabs(RefNucMap::angle(ca, frags[j].end_pos) - psi2);
            if (psi == RefNucMap::NULL_ANGLE) continue;
            psi = 1.0f - std::fmin(psi, 2.0f - psi);

            // Compute total score
            pss[k] = 0.5f * phi + 0.5f * psi;
        }
        return pss;
    }

    float phase_shift_score_wrt_reference(uint8_t chr_id, unsigned long x) const {
        float phi = this->angle(chr_id, x);
        if (phi == RefNucMap::NULL_ANGLE) {
            return RefNucMap::NULL_SCORE;
        } else {
            if (!distance) {
                return std::fmin(phi, 2.0f - phi);
            } else {
                return phi;
            }
        }
    }

    float alt_phase_shift_score_wrt_reference(uint8_t chr_id, unsigned long x, uint16_t frag_size) const {
        float phi = this->angle(chr_id, x);
        float psi = this->angle(chr_id, x + 1);
        if ((phi == RefNucMap::NULL_ANGLE) || (psi == RefNucMap::NULL_ANGLE)) {
            return RefNucMap::NULL_SCORE;
        } else {
            float nuc_distance = 2.0f / std::fabs(phi - psi);
            if (nuc_distance < 120.0f) nuc_distance = 120.0f;
            else if (nuc_distance > 250.0f) nuc_distance = 200.0f;
            return ((float) frag_size > nuc_distance);
        }
    }
};

} // namespace modalib

#endif // MODALIB_NUCLEOSOMES_H
