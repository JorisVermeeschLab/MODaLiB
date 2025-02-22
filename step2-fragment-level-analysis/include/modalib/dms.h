//  sample.h
//
//  Copyright 2024 Antoine Passemiers <antoine.passemiers@gmail.com>
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

#ifndef MODALIB_DMS_H
#define MODALIB_DMS_H

#include "sample.h"
#include "methylation.h"
#include "utils.h"
#include "array_utils.h"
#include "ordered_queue_map.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>

namespace modalib {

std::pair<float, size_t> compute_pearson_correlation(
        const std::vector<OrderedQueueMap<genome_size_t, CytosineContext>>& queues,
        genome_size_t c_coord,
        const std::vector<float>& y) {

    // Compute mean of the response variable and methylation mean across samples
    float c_mean = 0.0f;
    float y_mean = 0.0f;
    size_t n_samples_involved = 0;
    for (size_t i = 0; i < queues.size(); i++) {
        const OrderedQueueMap<genome_size_t, CytosineContext>& queue = queues.at(i);
        if (queue.contains(c_coord)) {
            const CytosineContext& c = queue.get(c_coord);
            c_mean += c.meth_ratio();
            y_mean += y.at(i);
            n_samples_involved++;
        }
    }
    if (n_samples_involved == 0) return {0.0f, 0};
    c_mean /= n_samples_involved;
    y_mean /= n_samples_involved;

    float var_y = 0.0f;
    float var_c = 0.0f;
    float cov = 0.0f;
    for (size_t i = 0; i < queues.size(); i++) {
        const OrderedQueueMap<genome_size_t, CytosineContext>& queue = queues.at(i);
        if (queue.contains(c_coord)) {
            const CytosineContext& c = queue.get(c_coord);
            float dev_c = c.meth_ratio() - c_mean;
            float dev_y = y.at(i) - y_mean;
            var_c += dev_c * dev_c;
            var_y += dev_y * dev_y;
            cov += dev_c * dev_y;
        }
    }
    float r = cov / std::sqrt(var_y * var_c);
    return {r, n_samples_involved};
}

void find_differentially_methylated_sites(
        const std::vector<std::string>& filepaths,
        const std::vector<float>& y,
        const std::string& out_filepath) {

    std::ofstream file;
    file.open(out_filepath);

    std::vector<OrderedQueueMap<genome_size_t, CytosineContext>> queues;
    std::vector<_Sample> samples;
    samples.reserve(filepaths.size());
    for (const std::string& filepath : filepaths) {
        samples.emplace_back(filepath);
        queues.emplace_back(50000);
    }
    size_t n = samples.size();

    // Store information about cytosines in array
    std::valarray<uint8_t> chr_ids(n);
    std::valarray<uint32_t> coords(n);
    std::valarray<genome_size_t> combined_coords(n);
    for (size_t i = 0; i < n; ++i) {
        const _Sample& sample = samples.at(i);
        const CytosineContext& c = sample.get_c();
        chr_ids[i] = (sample.eof()) ? 0 : c.chr_id;
        coords[i] = (sample.eof()) ? END_OF_SAMPLE : c.pos;
        combined_coords[i] = combine_coords(chr_ids[i], coords[i]);
    }

    while (true) {

        const int i = argmin(combined_coords);
        if (coords[i] == END_OF_SAMPLE) break;

        _Sample &sample = samples.at(i);

        const CytosineContext& c = sample.get_c();
        const uint8_t chr_id = c.chr_id;
        const uint32_t coord_i = c.pos;
        genome_size_t c_coords = combine_coords(chr_id, coord_i);

        // Compute Pearson correlation coefficient and store it
        OrderedQueueMap<genome_size_t, CytosineContext>& queue = queues.at(i);
        queue.add(c_coords, c);
        while (queue.size_exceeded() && (queue.size() > 0)) {
            genome_size_t pos = queue.front().first;
            const CytosineContext& c2 = queue.front().second;
            std::pair<float, size_t> res = compute_pearson_correlation(queues, pos, y);
            // std::cout << res.first << '\t' << res.second << std::endl;
            if ((res.second >= 50) && (std::fabs(res.first) > 0.4f)) {
                file << 1 * c2.chr_id << '\t' << c2.pos << '\t' << res.first << '\t' << res.second << '\t';
                if (c2.type == 1) file << "CpG";
                else if (c2.type == 2) file << "CHG";
                else if (c2.type == 3) file << "CHH";
                else file << "UNK";
                file << '\t';

                if (c2.type == 1) {
                    file << 'B';
                } else if ((c2.strand & 1) == 1) {
                    file << '+';
                } else {
                    file << '-';
                }

                file << '\n';
            }

            for (OrderedQueueMap<genome_size_t, CytosineContext>& q : queues) {
                if (q.contains(pos)) q.remove(pos);
            }
        }

        // Move to next fragment
        sample.next();

        // Update fragment statistics
        chr_ids[i] = (sample.eof()) ? 0 : c.chr_id;
        coords[i] = (sample.eof()) ? END_OF_SAMPLE : c.pos;
        combined_coords[i] = combine_coords(chr_ids[i], coords[i]);
    }

    file.close();
}

}  // namespace modalib

#endif //MODALIB_DMS_H
