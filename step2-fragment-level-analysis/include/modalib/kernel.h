//  kernel.h
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

#ifndef _KERNEL_H__
#define _KERNEL_H__

#include "modalib/counter.h"
#include "modalib/read_info.h"
#include "modalib/sample.h"
#include "modalib/osc.h"
#include "modalib/array_utils.h"
#include "bounds.h"
#include "utils.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <cstdint>
#include <functional>
#include <limits>
#include <math.h>
#include <unordered_map>
#include <utility>
#include <valarray>
#include <vector>

namespace modalib {

const long MAX_READ_DISTANCE = 150;
const long MAX_FRAGMENT_SIZE_DIFFERENCE = 1000;
const long MIN_FRAGMENT_SIZE = 20;


class KernelMatrix {
private:
    size_t n;
    std::valarray<long double> data;
    std::valarray<long double> counts;
public:
	KernelMatrix() : n(0), data(0), counts(0) {}
	KernelMatrix(size_t n): n(n), data(0.0, n * n), counts(0.0, n * n) {}
	KernelMatrix(const KernelMatrix& other): n(other.n), data(other.data), counts(other.counts) {}
	double operator()(int i, int j) noexcept {
		long double count = this->counts[i * this->n + j];
		if (count == 0) return 0;
		return this->data[i * this->n + j] / count;
	}
	inline void add_sym(float value, int i, int j, float weight) noexcept {
        this->counts[i * this->n + j] += weight;
        this->counts[j * this->n + i] += weight;
        long double vw = static_cast<long double>(value * weight);
        this->data[i * this->n + j] += vw;
        this->data[j * this->n + i] += vw;
	}
	void add_sym(
	        const std::valarray<float>& values, int i, const std::valarray<float>& weights, const std::vector<int>& idx
	        ) noexcept {
	    for (size_t k = 0; k < idx.size(); k++) {
	        const int j = idx[k];
	        this->counts[i * this->n + j] += weights[j];
	        this->counts[j * this->n + i] += weights[j];
            long double vw = static_cast<long double>(values[k] * weights[j]);
	        this->data[i * this->n + j] += vw;
	        this->data[j * this->n + i] += vw;
	    }
	}
	inline size_t n_rows() { return this->n; }
	inline size_t n_cols() { return this->n; }
};

void init_counter(const _Sample& _sample, Counter& counter) {
    _Sample sample(_sample);
    while (!sample.eof()) {
        const ReadInfo& read = sample.get();
        counter.inc_total_counts(read);
        sample.next();
    }
}

std::vector<Counter> init_counters(std::vector<_Sample>& samples) {
    const size_t n = samples.size();
    std::vector<Counter> counters(n);
    std::vector<std::thread> threads;
    for (size_t i = 0; i < n; i++) {
        threads.push_back(std::thread(init_counter, std::ref(samples.at(i)), std::ref(counters.at(i))));
    }
    for (std::thread& t : threads) t.join();
    return counters;
}

void compute_kernel_matrices_for_one_chromosome(
        int partition_id,
        const std::vector<_Sample>& _samples,
        std::unordered_map<uint8_t, std::vector<KernelMatrix>>& all_matrices,
        const Oscillator& osc,
        const std::vector<Counter>& _counters
        ) {

    size_t n = _samples.size();

    // Copy fragment counters -> one counter per sample and chromosome
    std::vector<Counter> counters;
    for (Counter _counter : _counters) {
        counters.push_back(Counter(_counter));
    }

    // Reset samples
    std::vector<_Sample> samples;
    for (size_t i = 0; i < n; ++i) {
        _Sample sample(_samples.at(i));
        sample.reset();

        // Only iterate over the fragments located on requested chromosome
        const std::pair<long, long> bounds = counters.at(i).get_bounds(partition_id);
        sample.set_lower_bound(bounds.first);
        sample.set_upper_bound(bounds.second);
        std::cout << "Bounds for sample " << i << ", partition " << partition_id << ": ";
        std::cout << bounds.first << " - " << bounds.second << std::endl;

        samples.push_back(sample);
    }

    // Indices of the samples whose current fragments are being compared to
    // current fragment of sample i.
    std::vector<int> idx(0);
    idx.reserve(n);

    std::valarray<ReadInfo> reads(n);
    std::valarray<uint8_t> chr_ids(n);
	std::valarray<long> coords(n);
	std::valarray<long> combined_coords(n);
	std::valarray<float> gc_ratios(n);
	for (size_t i = 0; i < n; ++i) {
	    _Sample& sample = samples.at(i);
	    const ReadInfo& read_i = sample.get();
	    reads[i] = read_i;
	    chr_ids[i] = (sample.eof()) ? 0 : read_i.start_chr_id;
		coords[i] = (sample.eof()) ? END_OF_SAMPLE : read_i.start_pos;
		combined_coords[i] = combine_coords(chr_ids[i], coords[i]);
		gc_ratios[i] = static_cast<float>((sample.eof()) ? 0 : read_i.get_gc_ratio());
	}

    unsigned int n_processed = 0;
	while (true) {

		const int i = argmin(combined_coords);
		if (coords[i] == END_OF_SAMPLE) break;

		_Sample& sample = samples.at(i);

        const ReadInfo& read_i = sample.get();
        const ReadInfo& previous_read_i = sample.get_previous();
        reads[i] = previous_read_i;
        const uint8_t chr_id = read_i.start_chr_id;

		if ((read_i.size >= MIN_FRAGMENT_SIZE) && (read_i.start_chr_id == read_i.end_chr_id)) {

            // Get kernel matrices for current chromosome
            std::vector<KernelMatrix>& matrices = all_matrices.at(chr_id);

            // Increment fragment counts for current sample i
            counters.at(i).add(read_i);

            // To compute k(x_i, x_i), don't compare current fragment with itself
            // but with previous fragment, to account for genome heterogeneity and noise
            chr_ids[i] = previous_read_i.start_chr_id;
            coords[i] = previous_read_i.start_pos;
            gc_ratios[i] = previous_read_i.get_gc_ratio();
            reads[i] = previous_read_i;

            // Get the indices of the samples whose current fragments are comparable with
            // current fragment of sample i
            idx.clear();
            const long coord_i = read_i.start_pos;
            std::valarray<float> w = read_i.get_gc_ratio() - gc_ratios;
            for (size_t j = 0; j < n; ++j) {

                // Alleviate sequencing biases by computing inter-fragment weights
                w[j] = 1.0f - std::fabs(w[j]);

                // The two fragments being compared should be sufficiently close
                w[j] *= (chr_id == chr_ids[j]);
                w[j] *= (std::labs(coord_i - coords[j]) <= MAX_READ_DISTANCE);

                // If the two fragments are comparable, then add sample j to the list
                // of valid indices
                if (w[j] > 0) idx.push_back(j);
            }

            // Update nucleosome positioning kernel matrix
            std::valarray<float> scores = osc.phase_shift_scores(chr_id, reads, idx, read_i);
            matrices.at(0).add_sym(scores, i, w, idx);

            const uint16_t l_i = read_i.size;
            const float m_i = read_i.get_meth_ratio();
            const float depth_i = static_cast<float>(read_i.n_sites());
            const float density_i = counters.at(i).get_density(read_i);
            const float mb_i = read_i.get_mutation_burden();
            for (size_t k = 0; k < idx.size(); ++k) {
                const int j = idx[k];

                // Update fragment length kernel matrix
                float score = ((l_i <= 150) == (reads[j].size <= 150));
                matrices.at(1).add_sym(score, i, j, w[j]);

                // Update methylation kernel matrix
                score = 1.0 - std::fabs(m_i - reads[j].get_meth_ratio());
                matrices.at(2).add_sym(score, i, j, w[j] * std::sqrt(depth_i * static_cast<float>(reads[j].n_sites())));

                // Update fragment counts kernel matrix
                const float density_j = counters.at(j).get_density(read_i);
                if (density_i * density_j > 0) {
                    const float den = std::fmax(density_i, density_j);
                    score = (den - std::fabs(density_i - density_j)) / den;
                    matrices.at(3).add_sym(score, i, j, 1.0f);
                }

                // Update "mutational" burden kernel matrix
                score = 1.0f - std::fabs(mb_i - reads[j].get_mutation_burden());
                matrices.at(4).add_sym(score, i, j, w[j] * std::sqrt(read_i.s_count * reads[j].s_count));
            }

        }

        n_processed++;

		// Move to next fragment
		sample.next();

		// Update fragment statistics
		chr_ids[i] = (sample.eof()) ? 0 : read_i.start_chr_id;
		coords[i] = (sample.eof()) ? END_OF_SAMPLE : read_i.start_pos;
		combined_coords[i] = combine_coords(chr_ids[i], coords[i]);
		gc_ratios[i] = static_cast<float>((sample.eof()) ? 0 : read_i.get_gc_ratio());
		reads[i] = read_i;
        
	}

	std::cout << "Finished processing partition " << partition_id << "." << std::endl;
    std::cout << "Number of fragments processed: " << n_processed << "." << std::endl;
}

std::unordered_map<uint8_t, std::vector<KernelMatrix>> compute_kernel_matrices(
        std::vector<_Sample>& samples,
        std::vector<bool>& chr_mask,
        std::unordered_map<uint8_t, std::vector<long>>& nuc_pos
        ) {

	size_t n = samples.size();

    // Initialize kernel matrices
    std::cout << "Initializing kernel matrices... ";
	std::unordered_map<uint8_t, std::vector<KernelMatrix>> matrices;
	for (uint8_t chr_id = 0; chr_id < chr_mask.size(); ++chr_id) {
	    if (chr_id < CHROMOSOME_BOUNDS.size()) {
            if (chr_mask.at(chr_id)) {
                matrices[chr_id] = std::vector<KernelMatrix>({
                    KernelMatrix(n),
                    KernelMatrix(n),
                    KernelMatrix(n),
                    KernelMatrix(n),
                    KernelMatrix(n)
                });
            }
        }
	}

    // Initialize fragment counters
    std::cout << "Computing total counts... ";
	std::vector<Counter> counters = init_counters(samples);
	std::cout << "done." << std::endl;

	// Initialize reference nucleosome map
	Oscillator osc(nuc_pos);

    // Compute kernel matrices
    std::vector<std::thread> threads;
	for (size_t partition_id = 0; partition_id < NUM_PARTITIONS; ++partition_id) {
        threads.push_back(std::thread(
            compute_kernel_matrices_for_one_chromosome,
            partition_id,
            std::ref(samples),
            std::ref(matrices),
            std::ref(osc),
            std::ref(counters)
        ));
	}
	for (std::thread& t : threads) t.join();

	return matrices;
}

} // namespace modalib

#endif // _KERNEL_H__
