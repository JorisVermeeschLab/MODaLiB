#ifndef MODALIB_FEATURES_H
#define MODALIB_FEATURES_H

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <valarray>
#include <vector>
#include <cstdlib>
#include <random>
#include <cmath>
#include <chrono>

#include "array_utils.h"
#include "fragment_info.h"
#include "nucleosomes.h"
#include "bounds.h"
#include "sample.h"
#include "region_list.h"
#include "utils.h"
#include "overlap.h"

namespace modalib {

const uint16_t MAX_FRAGMENT_SIZE = 1000;

static constexpr size_t IDX_NUCLEOSOME = 0;
static constexpr size_t IDX_FRAG_SIZE_INDEX = 1;
static constexpr size_t IDX_FRAG_LENGTH_COUNTS = 2;
static constexpr size_t IDX_METH = 3;
static constexpr size_t IDX_METH_SQUARED = 4;
static constexpr size_t IDX_CH_METH = 5;
static constexpr size_t IDX_READ_COUNTS = 6;
static constexpr size_t IDX_GC_CONTENT = 7;
static constexpr size_t IDX_CPG_COUNT = 8;
static constexpr size_t IDX_END_MOTIF_5P = 9;
static constexpr size_t IDX_END_MOTIF_3P = 10;
static constexpr size_t IDX_MAPPABILITY = 11;
static constexpr size_t IDX_NUC_SCORE_COUNTS = 12;
static constexpr size_t IDX_NUC_SCORE_DIST = 13;
static constexpr size_t IDX_FULLY_METH_FRAG_RATIO = 14;
static constexpr size_t IDX_CAFF = 15;
static constexpr size_t IDX_R1_METH = 16;
static constexpr size_t IDX_R2_METH = 17;
static constexpr size_t IDX_END_MOTIF_5P_171BP = 18;
static constexpr size_t IDX_LONG_FRAG_RATIO = 19;


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
    FeatureMatrix(const FeatureMatrix& other) = default;
    double operator()(size_t i, size_t j) const {
        long double count = this->counts[i * this->m + j];
        if (normalize) {
            if (count == 0) return 0.0;
            return static_cast<double>(this->data[i * this->m + j] / count);
        } else {
            return static_cast<double>(count);
        }
    }
    inline void add(float value, size_t i, size_t j, float weight) {
        this->counts[i * this->m + j] += weight;
        this->data[i * this->m + j] += value * weight;
    }
    inline void add_unweighted(float value, size_t i, size_t j, float weight) {
        this->counts[i * this->m + j] += weight;
        this->data[i * this->m + j] += value;
    }
    size_t n_rows() const { return this->n; }
    size_t n_cols() const { return this->m; }
};

std::vector<FeatureMatrix> extract_all_features(
        std::vector<_Sample>& samples,
        const RefNucMap& nuc_map,
        const _RegionList& _regions,
        double target_n_fragments,
        bool weighting
        ) {

    // Initialize list of regions
    _RegionList regions(_regions);
    regions.reset();

    // Initialize feature matrices
    size_t n = samples.size();
    size_t m = regions.get_num_regions();
    std::vector<FeatureMatrix> matrices({
        FeatureMatrix(n, m),  // Nucleosome positioning scores
        FeatureMatrix(n, m),  // Fragment size index
        FeatureMatrix(n, 500, false),  // Fragment length counts
        FeatureMatrix(n, m),  // CpG methylation (E[X])
        FeatureMatrix(n, m),  // CpG methylation (E[X^2])
        FeatureMatrix(n, m),  // CH methylation (E[X])
        FeatureMatrix(n, m, false),  // Frag counts
        FeatureMatrix(n, m),  // GC-content (E[X])
        FeatureMatrix(n, m, false),  // CpG counts
        FeatureMatrix(n, 256, false),  // 5' end-motifs
        FeatureMatrix(n, 256, false),  // 3' end-motifs
        FeatureMatrix(n, m),  // Mappability
        FeatureMatrix(n, m, false),  // Nuc score counts
        FeatureMatrix(n, 301, false),  // Nuc score histogram
        FeatureMatrix(n, m),  // Ratio of fully methylated fragments
        FeatureMatrix(n, m, false),  // Chromosomal aneuploidy of featured fragments
        FeatureMatrix(n, 500),  // R1 methylation
        FeatureMatrix(n, 500),  // R2 methylation
        FeatureMatrix(n, 256, false),  // 5' end-motifs (<= 171 bp)
        FeatureMatrix(n, m)  // Ratios of long fragments
    });

    // Add priors
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            matrices.at(IDX_NUCLEOSOME).add(0.500f, i, j, 0.005f);

            // Average methylation in human genome
            matrices.at(IDX_METH).add(0.787f, i, j, 0.001f);
            matrices.at(IDX_METH_SQUARED).add(0.787f * 0.787f, i, j, 0.005f);
            matrices.at(IDX_CH_METH).add(0.0f, i, j, 0.001f);

            // Average GC-content of human genome
            matrices.at(IDX_GC_CONTENT).add(0.41f, i, j, 0.005f);

            // Mappability
            matrices.at(IDX_MAPPABILITY).add(0.0f, i, j, 0.001f);
        }
    }

    for (size_t i = 0; i < n; ++i) {

        _Sample& sample = samples.at(i);

        // Initialize everything for in silico downsampling
        double downsampling_factor = 1.0;
        if (sample.get_length() > target_n_fragments) {
            downsampling_factor = static_cast<double>(target_n_fragments) / sample.get_length();
        }
        std::mt19937_64 rng;
        uint64_t seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::seed_seq ss{uint32_t(seed & 0xffffffff), uint32_t(seed >> 32)};
        rng.seed(ss);
        std::uniform_real_distribution<double> uniform_dist(0, 1);

        // Read the list of (marker) regions from the start
        regions.reset();

        while (!sample.eof()) {

            const FragmentInfo& frag_i = sample.get();
            const uint8_t chr_id = frag_i.chr_id;
            // const uint16_t l_i = frag_i.size;
            const uint16_t l_i = frag_i.insert_size;

            if (uniform_dist(rng) > downsampling_factor) {
                sample.next();  // Move to next fragment
                continue;
            }

            if (frag_i.mapping_quality < 45) {  // TODO: remove
                sample.next();  // Move to next fragment
                continue;
            }

            float gc_weight = 1.0f;
            if (weighting) gc_weight = frag_i.weight;

            // Filter out fragments that are not in the region list provided by the user.
            // In order to speed up computations by reducing the search time of this method,
            // the input file is assumed to be "almost" sorted by chromosomal coordinates.
            Overlap ov = regions.search_overlap(frag_i);
            while (ov.fraction > 0) {

                // Get region index
                const size_t j = ov.region_id;

                // Update statistics about end-motifs
                if (l_i <= MAX_FRAGMENT_SIZE) {
                    if (ov.left_end_inside()) {
                        matrices.at(IDX_END_MOTIF_5P).add(1.0f, i, frag_i.end_motif_5p, gc_weight);
                        if (frag_i.strand % 2 == 0) {
                            matrices.at(IDX_END_MOTIF_5P_171BP).add(1.0f, i, frag_i.end_motif_5p, gc_weight);
                        }
                    }
                    if (ov.right_end_inside()) {
                        matrices.at(IDX_END_MOTIF_3P).add(1.0f, i, frag_i.end_motif_3p, gc_weight);
                    }
                }

                // Update nucleosome features
                if (ov.left_end_inside()) {
                    float score1 = nuc_map.phase_shift_score_wrt_reference(chr_id, static_cast<unsigned long>(frag_i.start_pos));
                    if (score1 != RefNucMap::NULL_SCORE) {
                        matrices.at(IDX_NUCLEOSOME).add(score1, i, j, 0.5f * gc_weight);
                        matrices.at(IDX_NUC_SCORE_COUNTS).add(1.0f, i, j, gc_weight);
                        size_t nuc_i = static_cast<size_t>(std::roundf(score1 * 300.0f));
                        matrices.at(IDX_NUC_SCORE_DIST).add(1.0f, i, nuc_i, gc_weight);
                    }
                }
                if (ov.right_end_inside()) {
                    float score2 = nuc_map.phase_shift_score_wrt_reference(chr_id, static_cast<unsigned long>(frag_i.end_pos));
                    if (score2 != RefNucMap::NULL_SCORE) {
                        matrices.at(IDX_NUCLEOSOME).add(score2, i, j, 0.5f * gc_weight);
                        matrices.at(IDX_NUC_SCORE_COUNTS).add(1.0f, i, j, gc_weight);
                        size_t nuc_i = static_cast<size_t>(std::roundf(score2 * 300.0f));
                        matrices.at(IDX_NUC_SCORE_DIST).add(1.0f, i, nuc_i, gc_weight);
                    }
                }

                // Update proportion of long fragments
                if (((100 <= l_i) && (l_i <= 166)) || ((169 <= l_i) && (l_i <= 240))) {
                    float length_score = (l_i >= 169) ? 1.0f : 0.0f;
                    matrices.at(IDX_FRAG_SIZE_INDEX).add(length_score, i, j, ov.fraction * gc_weight);
                }
                float length_score2 = (l_i >= 167) ? 1.0f : 0.0f;
                    matrices.at(IDX_LONG_FRAG_RATIO).add(length_score2, i, j, ov.fraction * gc_weight);

                // Update fragment length counts
                if ((l_i <= MAX_FRAGMENT_SIZE) && (frag_i.mapping_quality >= 60)) {  // TODO
                    int hist_max_index = matrices.at(IDX_FRAG_LENGTH_COUNTS).n_cols() - 1;
                    int k = (l_i < hist_max_index) ? l_i : hist_max_index;
                    matrices.at(IDX_FRAG_LENGTH_COUNTS).add(1.0f, i, k, ov.fraction * gc_weight);
                }

                // Update genome-wide methylation features
                const float m_i = frag_i.get_cpg_meth_ratio();
                matrices.at(IDX_METH).add(m_i, i, j, frag_i.n_cpg_sites() * ov.fraction * gc_weight);
                matrices.at(IDX_METH_SQUARED).add(m_i * m_i, i, j, frag_i.n_cpg_sites() * ov.fraction * gc_weight);
                matrices.at(IDX_CH_METH).add(frag_i.get_ch_meth_ratio(), i, j, frag_i.n_ch_sites() * ov.fraction * gc_weight);
                matrices.at(IDX_FULLY_METH_FRAG_RATIO).add((m_i == 1.0f), i, j, frag_i.n_cpg_sites() * ov.fraction * gc_weight);
                if (l_i <= MAX_FRAGMENT_SIZE) {
                    uint16_t hist_max_index = matrices.at(IDX_FRAG_LENGTH_COUNTS).n_cols() - 1;
                    int k = (l_i < hist_max_index) ? l_i : hist_max_index;
                    matrices.at(IDX_R1_METH).add(frag_i.get_r1_meth_ratio(), i, k, frag_i.n_cpg_sites() * ov.fraction * gc_weight);
                    matrices.at(IDX_R2_METH).add(frag_i.get_r2_meth_ratio(), i, k, frag_i.n_cpg_sites() * ov.fraction * gc_weight);
                }

                // Update fragment counts
                matrices.at(IDX_READ_COUNTS).add(1.0f, i, j, ov.fraction * gc_weight);
                if ((l_i < 151) || (l_i > 220)) {
                    matrices.at(IDX_CAFF).add(1.0f, i, j, ov.fraction * gc_weight);
                }

                // Update extra info (GC-content, number of CpG sites, etc.)
                const float gc_i = frag_i.get_gc_ratio();
                matrices.at(IDX_GC_CONTENT).add(gc_i, i, j, frag_i.s_count * ov.fraction * gc_weight);
                matrices.at(IDX_CPG_COUNT).add(1.0f, i, j, frag_i.n_cpg_sites() * ov.fraction * gc_weight);

                // Update average mappability score
                matrices.at(IDX_MAPPABILITY).add(frag_i.mapping_quality, i, j, ov.fraction * gc_weight);

                // If the fragment is not overlapping with the next region,
                // then move to the next fragment.
                if (!ov.may_overlap_next_region) {
                    break;
                } else {
                    ov = regions.overlap(frag_i, j + 1);
                }
            }

            // Move to next fragment
            sample.next();
        }
    }

    return matrices;
}

} // namespace modalib

#endif // MODALIB_FEATURES_H
