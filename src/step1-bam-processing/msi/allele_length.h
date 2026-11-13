#ifndef MODALIB_MSI_ALLELE_LENGTH_H
#define MODALIB_MSI_ALLELE_LENGTH_H

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "../read_info.h"
#include "vntr.h"

// A locus is located on the reference; observed_length is measured in the
// sample's aligned query bases. Coordinates are 1-based inclusive here.
struct MicrosatelliteAllele {
    uint32_t start = 0;
    uint32_t end = 0;
    uint32_t reference_length = 0;
    uint32_t observed_length = 0;
    std::string motif;

    // Fractional copy counts are meaningful as length / motif period;
    // they do not, by themselves, demonstrate an uninterrupted motif sequence.
    double reference_repeat_count() const {
        return motif.empty() ? 0.0 :
            static_cast<double>(reference_length) / static_cast<double>(motif.size());
    }
    double observed_repeat_count() const {
        return motif.empty() ? 0.0 :
            static_cast<double>(observed_length) / static_cast<double>(motif.size());
    }
};

inline bool msi_is_canonical(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

// Require aligned query bases on both sides of a reference repeat.
// Do not count an insertion as an anchor. Anchor quality matters; converted
// bases themselves are allowed, since the sequencing chemistry is EM-seq.
inline bool msi_has_anchors(const ReadInfo& read, uint32_t first,
                            uint32_t last, uint32_t flank,
                            uint8_t min_phred = 5) {
    if (first == 0 || last < first || first <= flank ||
        static_cast<uint64_t>(last) + flank > UINT32_MAX) return false;
    if (flank == 0) return true;
    auto anchor = [&](uint32_t position) {
        auto it = std::lower_bound(read.bases.begin(), read.bases.end(), position,
           [](const std::pair<uint32_t, BaseInfo>& p, uint32_t pos) {
               return p.first < pos;
           });
        for (; it != read.bases.end() && it->first == position; ++it) {
            const BaseInfo& b = it->second;
            if (!b.inserted && !b.deleted && !b.soft_clipped &&
                msi_is_canonical(b.c) && b.phred_score >= min_phred)
                return true;
        }
        return false;
    };
    for (uint32_t d = 1; d <= flank; ++d) {
        if (!anchor(first - d) || !anchor(last + d)) return false;
    }
    return true;
}

// Count query-consuming bases aligned inside the REFERENCE locus.
// A sample insertion following reference position p shares p's coordinate in
// process_read(), so it counts when first <= p < last. Deletions count zero.
// Reject skipped/uncovered reference positions and ambiguous structural spans.
inline bool msi_observed_span(const ReadInfo& read, uint32_t first,
                              uint32_t last, uint32_t& observed_length) {
    if (!first || last < first) return false;
    uint64_t count = 0;
    uint64_t expected = first;
    for (const auto& entry : read.bases) {
        const uint32_t pos = entry.first;
        const BaseInfo& b = entry.second;
        if (pos < first) continue;
        if (pos > last) break;
        if (b.soft_clipped) return false;
        if (b.inserted) {
            if (pos < last) ++count;  // Not an insertion after the right boundary
            continue;
        }
        if (pos != expected) return false;
        ++expected;
        if (!b.deleted) ++count;
    }
    if (expected != static_cast<uint64_t>(last) + 1 ||
        count > UINT32_MAX) return false;
    observed_length = static_cast<uint32_t>(count);
    return true;
}

// Find exact tandem repeats in the reference bases represented by this read,
// then estimate the observed allele span using the actual aligned read bases.
// A return value is only a putative repeat-length allele: it does not prove
// that the converted sample sequence still consists of exact motif copies.
inline std::vector<MicrosatelliteAllele> msi_call_alleles(
    const ReadInfo& read, uint32_t flank = 3, uint8_t min_phred = 5) {

    std::vector<MicrosatelliteAllele> calls;
    if (read.bases.size() != read.ref_bases.size()) return calls;

    std::vector<std::pair<uint32_t, BaseInfo>> reference_bases;
    reference_bases.reserve(read.bases.size());
    for (size_t i = 0; i < read.bases.size(); ++i) {
        const auto& observed = read.bases[i];
        const BaseInfo& ref = read.ref_bases[i].second;
        if (observed.second.inserted || observed.second.soft_clipped ||
            ref.deleted || ref.soft_clipped) continue;
        // Retain reference bases even where the sample has a deletion.
        reference_bases.emplace_back(observed.first, ref);
    }

    std::vector<TandemRepeat> reference_repeats;
    find_microsatellites(reference_bases, reference_repeats);
    for (const auto& tr : reference_repeats) {
        if (!tr.is_valid || tr.start < 1 || tr.end < tr.start) continue;
        if (!msi_has_anchors(read, tr.start, tr.end, flank, min_phred)) continue;
        uint32_t observed = 0;
        if (!msi_observed_span(read, tr.start, tr.end, observed)) continue;
        if (!observed) continue;
        const uint32_t reference_length = tr.end - tr.start + 1;
        calls.push_back({tr.start, tr.end, reference_length, observed, tr.motif});
    }
    return calls;
}

#endif // MODALIB_MSI_ALLELE_LENGTH_H
