//  microsatellites.h
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

#ifndef MODALIB_MICROSATELLITES_H
#define MODALIB_MICROSATELLITES_H

#include <cstdlib>
#include <map>
#include <utility>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <cmath>
#include <set>
#include <ctime>
#include <string>
#include <sstream>
#include <thread>
#include <unordered_set>

#include "features.h"
#include "utils.h"

namespace modalib {

class MSS {
private:
    std::string filepath;
    std::map<genome_size_t, std::unordered_map<int, size_t>> data;
    std::map<genome_size_t, size_t> counts;
    std::map<genome_size_t, std::string> motifs;
public:
    MSS() = default;
    explicit MSS(std::string filepath) : filepath(std::move(filepath)) {}
    ~MSS() = default;
    inline bool contains_ms_at(genome_size_t combined_coord);
    inline bool contains_ms_at(genome_size_t combined_coord, int indels);
    void load();
    void close();
    void extract_features(int i, FeatureMatrix& features, FeatureMatrix& counts);
    void extract_markers(int i, FeatureMatrix& features, FeatureMatrix& counts, const std::unordered_map<genome_size_t, size_t>& markers_pos);
};


void MSS::load() {
    data.clear();
    counts.clear();
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filepath);
    } else {
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;
            std::string motif;
            std::size_t i = 0;
            uint8_t chr_id = 0;
            uint32_t start = 0;
            int indels = 0;
            while (std::getline(ss, token, ',')) {
                switch(i) {
                    case 0:
                        chr_id = std::stoi(token);
                        break;
                    case 1:
                        start = std::stoi(token);
                        break;
                    case 3:
                        motif = token;
                        break;
                    case 4:
                        indels = std::stoi(token);
                        break;
                    default:
                        break;
                }
                i++;
            }
            if ((chr_id > 0) && (start > 0)) {
                genome_size_t combined_coords = combine_coords(chr_id, start);
                if (data.find(combined_coords) == data.end()) {
                    data[combined_coords] = std::unordered_map<int, size_t>();
                    counts[combined_coords] = 0;
                }
                if (data[combined_coords].find(combined_coords) == data[combined_coords].end()) {
                    data[combined_coords][indels] = 0;
                }
                data[combined_coords][indels]++;
                counts[combined_coords]++;
                motifs[combined_coords] = motif;
            }
        }
    }
}


inline bool MSS::contains_ms_at(genome_size_t combined_coords) {
    return this->data.find(combined_coords) != this->data.end();
}


inline bool MSS::contains_ms_at(genome_size_t combined_coords, int indels) {
    return this->data.at(combined_coords).find(indels) != this->data.at(combined_coords).end();
}


void MSS::close() {
    data.clear();
    counts.clear();
    motifs.clear();
}


inline size_t get_feature_index_from_motif(size_t n, size_t max_n_alleles) {
    return (n >= max_n_alleles) ? max_n_alleles - 1 : n;
}


void MSS::extract_features(int i, FeatureMatrix& features, FeatureMatrix& c) {

    size_t max_n_alleles = features.n_cols();

    // Read data from files
    this->load();

    for (const auto& it : this->data) {

        genome_size_t pos = it.first;

        // Get counts from current microsatellite
        size_t n = this->counts.at(pos);

        // Get feature index from motif
        size_t j = get_feature_index_from_motif(n, max_n_alleles);

        c.add(1.0f, i, j, 1.0f);

        float p;
        // if (n <= 2) {
        if (n <= 1) {  // TODO
            p = 0.0f;
        } else {

            n = 0;
            size_t m1 = 0;
            size_t m2 = 0;
            for (const auto &it2 : this->data.at(pos)) {
                if (it2.second > m1) {
                    m2 = m1;
                    m1 = it2.second;
                } else if (it2.second > m2) {
                    m2 = it2.second;
                }
                n += it2.second;
            }
            p = static_cast<float>(n - m1 - m2) / static_cast<float>(n);
        }

        features.add(p, i, j, 1.0f);
    }

    // Remove data from memory
    this->close();
}


void MSS::extract_markers(int i, FeatureMatrix& features, FeatureMatrix& c, const std::unordered_map<genome_size_t, size_t>& markers_pos) {

    // Read data from files
    this->load();

    for (const auto& it : this->data) {

        genome_size_t pos = it.first;

        // If microsatellite is not in marker list, skip it.
        if (markers_pos.find(pos) == markers_pos.end()) continue;

        // Get counts from current microsatellite
        size_t n = this->counts.at(pos);

        // Get feature index from motif
        size_t j = markers_pos.at(pos);

        if (j >= c.n_cols()) continue;

        c.add(n, i, j, 1.0f);

        float p;
        if (n <= 2) {
            p = 0.0f;
        } else {

            n = 0;
            size_t m1 = 0;
            size_t m2 = 0;
            for (const auto &it2 : this->data.at(pos)) {
                if (it2.second > m1) {
                    m2 = m1;
                    m1 = it2.second;
                } else if (it2.second > m2) {
                    m2 = it2.second;
                }
                n += it2.second;
            }
            p = static_cast<float>(n - m1) / static_cast<float>(n);
        }

        features.add(p, i, j, 1.0f);
    }

    // Remove data from memory
    this->close();
}


} // namespace modalib

#endif //MODALIB_MICROSATELLITES_H
