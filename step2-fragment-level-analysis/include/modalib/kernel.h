//  kernel.h
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

#ifndef MODALIB_KERNEL_H
#define MODALIB_KERNEL_H

#include "counter.h"
#include "fragment_info.h"
#include "sample.h"
#include "array_utils.h"
#include "bounds.h"
#include "features.h"
#include "nucleosomes.h"
#include "utils.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <cstdint>
#include <functional>
#include <limits>
#include <cmath>
#include <unordered_map>
#include <utility>
#include <valarray>
#include <vector>

namespace modalib {


// Number of kernel functions
const int NUM_KERNELS = 5;

// Maximum distance allowed between the start positions of two fragments.
// Above that threshold, the fragment pair will be ignored and won't contribution
// to the estimation of any of the kernel functions.
const long MAX_FRAGMENT_DISTANCE = 2000;

// Minimum fragment size.
// Below that threshold, the current fragment won't contribute to the estimation
// of any of the kernel functions.
const unsigned long MIN_FRAGMENT_SIZE = 0;


class KernelMatrix {
private:

    // Number of rows (= number of columns)
    size_t n;

    // Flattened symmetric matrix of shape n*n representing the unnormalized kernel matrix.
    // Kernel between samples S_i and S_j is computed as:
    // k(S_i, S_j) = `data[i*n+j]` / `weights[i*n+j]`
    std::valarray<long double> data;

    // Flattened symmetric matrix of shape n*n representing the normalization factors
    // for the kernel matrix.
    std::valarray<long double> weights;

public:

    // Default constructor (kernel matrix of size 0*0)
	KernelMatrix() : n(0), data(0), weights(0) {}

    // Constructor (kernel matrix of size n*n)
	KernelMatrix(size_t n): n(n), data(0.0, (n * (n + 1)) / 2), weights(0.0, (n * (n + 1)) / 2) {}

    // Copy constructor
	KernelMatrix(const KernelMatrix& other): n(other.n), data(other.data), weights(other.weights) {}

    /**
     * Converts a couple of flat indices to a flat index.
     * 
     * This function accounts for the symmetric structure of the matrix.
     *
     * @param i Row index
     * @param j Column index
     * @return Flat index.
     */
    [[nodiscard]] inline int convert_index(int i, int j) const {
        if (i > j) {
            int temp = i;
            i = j;
            j = temp;
        }
        return i * this->n + j - (i * (i + 1)) / 2;
    }

    /**
     * Access element (i, j) of the kernel matrix.
     *
     * @param i Row index
     * @param j Column index
     * @return Element (i, j) of the kernel matrix
     */
	inline double operator()(int i, int j) {
        const int k = this->convert_index(i, j);
		const long double weight = this->weights[k];
		if (weight == 0.0) return 0.0;
		return this->data[k] / weight;
	}

    /**
     * Add value to element (i, j) of the kernel matrix.
     *
     * Sample-level kernel matrix is defined as a weighted linear combination
     * of fragment-level kernels: \f[
     *     k(S_a, S_b) = \frac{\sum_{x_1 \in S_a}\sum_{x_2 \in S_a} w(x_1, x_2) \tilde{k}(x_1, x_2)}
     *         {\sum_{x_1 \in S_a}\sum_{x_2 \in S_a} w(x_1, x_2)}
     * \f].
     * Attribute `data` contains the numerator of \f$k(S_a, S_b)\f$ and attribute `weights` contains
     * the denominator of \f$k(S_a, S_b)\f$. Therefore, numerator and denominators should be update
     * separately, hence the need for an extra `weight` method parameter.
     *
     * @param value Value to multiply with `weight` and add (inplace) to the numerator of the kernel function.
     * @param i Row index
     * @param j Column index
     * @param weight Value to add (inplace) to the denominator of the kernel function.
     */
	inline void add_sym(float value, int i, int j, float weight) {
        const int k = this->convert_index(i, j);
        this->weights[k] += weight;
        this->data[k] += static_cast<long double>(value * weight);
	}

    // Number of rows (= number of columns)
	inline size_t n_rows() const { return this->n; }

    // Number of columns (= number of rows)
	inline size_t n_cols() const { return this->n; }
};


/**
 * Count the fragments in one sample.
 *
 * @param sample Sample to be processed
 * @param counter Counter where to store the counts per GC-stratum (unique to current sample).
 * @param counts Per-bin counts, shared across all samples.
 */
void init_counter(const _Sample& _sample, Counter& counter, std::unordered_map<uint8_t, FeatureMatrix>& counts) {
    _Sample sample(_sample);
    while (!sample.eof()) {
        const FragmentInfo& frag = sample.get();

        // Increment per-bin counts
        const uint8_t chr_id = frag.start_chr_id;
        const long coord = frag.start_pos;
        const size_t j = coord / BIN_SIZE;
        if (j < counts.at(chr_id).n_cols()) {
            // `counts` matrix is of shape (1, p), as the per-bin counts
            // are shared across all samples. For this reason, we increment the 1 first row only.
            counts.at(chr_id).add(1.0f, 0, j, 1.0f);
        }

        // Increment per-GC-stratum counts
        counter.inc_total_counts(frag);

        sample.next();
    }
}


/**
 * Count the fragments in each sample.
 *
 * Samples are processed in parallel using multithreading.
 *
 * @param sample Samples to be processed
 * @param counts Per-bin counts
 * @return Vector of counters containing the fragment counts.
 *     the vector has the same size as `samples`.
 */
std::vector<Counter> init_counters(std::vector<_Sample>& samples, std::unordered_map<uint8_t, FeatureMatrix>& counts) {
    const size_t n = samples.size();
    std::vector<Counter> counters(n);
    std::vector<std::thread> threads;
    for (size_t i = 0; i < n; i++) {
        threads.push_back(std::thread(
            init_counter, std::ref(samples.at(i)), std::ref(counters.at(i)), std::ref(counts)
        ));
    }
    for (std::thread& t : threads) t.join();
    return counters;
}


/**
 * Compute the kernel matrices for one partition.
 *
 * Each partition is processed by a different thread,
 * allowing the parallel processing of multiple (non-overlapping) portions of the genome.
 *
 * @param partition_id 0-based index of current partition.
 * @param _samples Samples for which to compute the kernel matrices.
 *     Assuming there are `n` samples, each kernel matrix will be of shape `n × n`.
 * @param all_matrices Kernel matrices to be computed.
 * TODO
 */
void compute_kernel_matrices_for_one_chromosome(
        int partition_id,
        const std::vector<_Sample>& _samples,
        std::unordered_map<uint8_t, std::vector<KernelMatrix>>& all_matrices,
        const RefNucMap& nuc_map,
        const std::vector<Counter>& _counters,
        const std::unordered_map<uint8_t, FeatureMatrix>& counts
        ) {

    size_t n = _samples.size();

    // Copy fragment counters -> one counter per sample and chromosome
    std::vector<Counter> counters;
    for (const Counter& _counter : _counters) {
        counters.emplace_back(_counter);
    }

    // Reset samples
    std::vector<_Sample> samples;
    for (size_t i = 0; i < n; ++i) {
        _Sample sample(_samples.at(i));

        // Only iterate over the fragments located on the partition associated
        // with current thread
        const std::pair<long, long> bounds = counters.at(i).get_bounds(partition_id);
        if (bounds.first == bounds.second) return;
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

    // Store information about fragments in arrays
    std::valarray<FragmentInfo> frags(n);
    std::valarray<uint8_t> chr_ids(n);
	std::valarray<long> coords(n);
    std::valarray<uint16_t> lengths(n);
	std::valarray<genome_size_t> combined_coords(n);
	std::valarray<float> gc_ratios(n);
	for (size_t i = 0; i < n; ++i) {
	    _Sample& sample = samples.at(i);
	    const FragmentInfo& frag_i = sample.get();
	    frags[i] = frag_i;
	    chr_ids[i] = (sample.eof()) ? 0 : frag_i.start_chr_id;
		coords[i] = (sample.eof()) ? END_OF_SAMPLE : frag_i.start_pos;
		combined_coords[i] = combine_coords(chr_ids[i], coords[i]);
        lengths[i] = (sample.eof()) ? 0 : frag_i.size;
        gc_ratios[i] = frag_i.get_gc_ratio();
	}

    unsigned int n_processed = 0;
	while (true) {

		const int i = argmin(combined_coords);
		if (coords[i] == END_OF_SAMPLE) break;

		_Sample& sample = samples.at(i);

        const FragmentInfo& frag_i = sample.get();
        const uint8_t chr_id = frag_i.start_chr_id;
        const long coord_i = frag_i.start_pos;

        // Determine kernel bin
        const unsigned int kernel_bin_id = coord_i / KERNEL_BIN_SIZE;

        // Increment fragment counts for current sample i
        counters.at(i).add(frag_i);

        // Compute bin weight
        float bin_weight = 0.0f;
        const int bin_id = coord_i / BIN_SIZE;
        if ((bin_id < counts.at(chr_id).n_cols()) && (counts.at(chr_id)(0, bin_id)) > 0) {
            bin_weight = 1.0 / counts.at(chr_id)(0, bin_id);
        }

        // Count chromosome mismatches (proxy for translocations)
        if (frag_i.start_chr_id != frag_i.end_chr_id) {
            // TODO
        }

        // Get kernel matrices for current chromosome
        std::vector<KernelMatrix>& matrices = all_matrices.at(chr_id);

		if ((frag_i.size >= MIN_FRAGMENT_SIZE) && (frag_i.start_chr_id == frag_i.end_chr_id) && (NUM_KERNELS * kernel_bin_id < matrices.size())) {

            // To compute k(x_i, x_i), don't compare current fragment with itself
            // but with previous fragment, to account for genome heterogeneity and noise
            const FragmentInfo& previous_frag_i = sample.get_previous();
            frags[i] = previous_frag_i;
            chr_ids[i] = previous_frag_i.start_chr_id;
            coords[i] = previous_frag_i.start_pos;
            combined_coords[i] = combine_coords(chr_ids[i], coords[i]);
            lengths[i] = previous_frag_i.size;
            gc_ratios[i] = previous_frag_i.get_gc_ratio();

            // Get the indices of the samples whose current fragments are comparable with
            // current fragment of sample i
            idx.clear();
            std::valarray<float> w = frag_i.get_gc_ratio() - gc_ratios;
            for (size_t j = 0; j < n; ++j) {

                // Alleviate sequencing biases by computing inter-fragment weights
                w[j] = 1.0f;
                // w[j] = 1.0f - std::fabs(w[j]);
                // w[j] *= w[j];

                // The two fragments being compared should be sufficiently close
                w[j] *= static_cast<float>(chr_id == chr_ids[j]);
                w[j] *= static_cast<float>(std::labs(coord_i - coords[j]) <= MAX_FRAGMENT_DISTANCE);

                // Account for the coverage -> downweight high-coverage bins and vice versa
                // w[j] *= bin_weight;

                // Give more importance to fragment pairs with high overlap
                w[j] *= 1.0f + 2.0f * frag_i.overlap(frags[j]);

                // If the two fragments are comparable, then add sample j to the list
                // of valid indices
                if (w[j] > 0) idx.push_back(j);

            }

            // Compute phase shift scores (needed to update the nucleosome positioning kernel matrix).
            // `scores` has the same size as `idx`.
            std::valarray<float> scores = nuc_map.phase_shift_scores(chr_id, frags, idx, frag_i);

            const uint16_t l_i = frag_i.size;
            const float m_i = frag_i.get_meth_ratio();
            const auto depth_i = static_cast<float>(frag_i.n_sites());
            const float density_i = counters.at(i).get_density(frag_i);
            const float mb_i = frag_i.get_mutation_burden();
            for (size_t k = 0; k < idx.size(); ++k) {
                const int j = idx[k];

                // Update nucleosome positioning kernel matrix
                matrices.at(NUM_KERNELS * kernel_bin_id).add_sym(scores[k], i, j, w[j]);

                // Update fragment length kernel matrix.
                // k(x_i, x_j) = 1 when both fragments are on the same side w.r.t.
                // the expected modal value (167 bp), and 0 otherwise.
                // In other words, k(x_i, x_j) = 1 when x_i and x_j are either both
                // "long" or both "short" fragments.
                // Weighting function is the default one.
                // float max_score = std::fmax(l_i, frags[j].size);
                // float score = (max_score - std::fabs(l_i - frags[j].size)) / max_score;
                float score = ((l_i <= 167) == (frags[j].size <= 167));  // TODO
                matrices.at(NUM_KERNELS * kernel_bin_id + 1).add_sym(score, i, j, w[j]);

                // Update methylation kernel matrix.
                score = 1.0f - std::fabs(m_i - frags[j].get_meth_ratio());
                matrices.at(NUM_KERNELS * kernel_bin_id + 2).add_sym(score, i, j, w[j] * std::sqrt(depth_i * static_cast<float>(frags[j].n_sites())));

                // Update fragment counts kernel matrix
                const float density_j = counters.at(j).get_density(frag_i);
                if (density_i * density_j > 0) {
                    const float den = std::fmax(density_i, density_j);
                    score = (den - std::fabs(density_i - density_j)) / den;
                    matrices.at(NUM_KERNELS * kernel_bin_id + 3).add_sym(score, i, j, bin_weight);
                }

                // Update "mutational" burden kernel matrix
                score = 1.0f - std::fabs(mb_i - frags[j].get_mutation_burden());
                matrices.at(NUM_KERNELS * kernel_bin_id + 4).add_sym(score, i, j, w[j] * std::sqrt(frag_i.s_count * frags[j].s_count));
            }
        }

        n_processed++;

		// Move to next fragment
		sample.next();

		// Update fragment statistics
		chr_ids[i] = (sample.eof()) ? 0 : frag_i.start_chr_id;
		coords[i] = (sample.eof()) ? END_OF_SAMPLE : frag_i.start_pos;
		combined_coords[i] = combine_coords(chr_ids[i], coords[i]);
        lengths[i] = (sample.eof()) ? 0 : frag_i.size;
		gc_ratios[i] = sample.get().get_gc_ratio();
		frags[i] = frag_i;
	}

	std::cout << "Finished processing partition " << partition_id << "." << std::endl;
    std::cout << "Number of fragments processed: " << n_processed << "." << std::endl;
}

std::unordered_map<uint8_t, std::vector<KernelMatrix>> compute_kernel_matrices(
        std::vector<_Sample>& samples,
        RefNucMap& nuc_map
        ) {

	size_t n = samples.size();

    // Initialize kernel matrices
    std::cout << "Initializing kernel matrices... ";
	std::unordered_map<uint8_t, std::vector<KernelMatrix>> matrices;
	for (uint8_t chr_id = 0; chr_id < CHROMOSOME_BOUNDS.size(); ++chr_id) {
        matrices[chr_id] = std::vector<KernelMatrix>();
        for (int k = 0; k < KERNEL_BOUNDS.at(chr_id).second * NUM_KERNELS; k++) {
            matrices[chr_id].push_back(KernelMatrix(n));
        }
	}

    // Initialize fragment counters
    std::unordered_map<uint8_t, FeatureMatrix> counts;
    for (uint8_t chr_id = 0; chr_id < CHROMOSOME_BOUNDS.size(); ++chr_id) {
        size_t m = CHROMOSOME_BOUNDS.at(chr_id).second;
        counts[chr_id] = FeatureMatrix(1, m, false);
    }
    std::cout << "Computing total counts... ";
	std::vector<Counter> counters = init_counters(samples, counts);
	std::cout << "done." << std::endl;

    // Compute kernel matrices
    std::vector<std::thread> threads;
	for (size_t partition_id = 0; partition_id < NUM_PARTITIONS; ++partition_id) {
        threads.push_back(std::thread(
            compute_kernel_matrices_for_one_chromosome,
            partition_id,
            std::ref(samples),
            std::ref(matrices),
            std::ref(nuc_map),
            std::ref(counters),
            std::ref(counts)
        ));
	}
	for (std::thread& t : threads) t.join();

	return matrices;
}

} // namespace modalib

#endif // MODALIB_KERNEL_H
