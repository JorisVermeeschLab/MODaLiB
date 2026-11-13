//  vntr.cpp -- coordinate-safe tandem repeat detection
//  Copyright 2023 Antoine Passemiers
//  SPDX-License-Identifier: GPL-2.0-or-later

#include "vntr.h"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace {

bool canonical(char c) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

// Report maximal, nonoverlapping exact tandem repeats. Coordinates are
// 0-based inclusive relative to `seq`.
void scan_sequence(const char* seq, size_t n, std::vector<TandemRepeat>& out) {
    if (!seq) return;
    size_t i = 0;
    while (i < n) {
        if (!canonical(seq[i])) { ++i; continue; }
        TandemRepeat best;
        size_t best_length = 0;
        for (int motif_len = MIN_MOTIF_LENGTH; motif_len <= MAX_MOTIF_LENGTH; ++motif_len) {
            const size_t k = static_cast<size_t>(motif_len);
            if (k > (n - i) / MIN_NUM_REPEATS) continue;
            bool motif_ok = true;
            for (size_t j = 0; j < k; ++j) {
                if (!canonical(seq[i+j])) { motif_ok = false; break; }
            }
            if (!motif_ok) continue;

            size_t copies = 1;
            while (copies < (n-i)/k) {
                const size_t next = i + copies*k;
                bool equal = true;
                for (size_t j = 0; j < k; ++j) {
                    if (!canonical(seq[next+j]) ||
                        std::toupper(static_cast<unsigned char>(seq[i+j])) !=
                            std::toupper(static_cast<unsigned char>(seq[next+j]))) {
                        equal = false;
                        break;
                    }
                }
                if (!equal) break;
                ++copies;
            }
            // Recognize terminal partial motif copies (e.g. ACT*4 + AC).
            // A terminal partial copy must agree with the beginning of the motif.
            size_t partial = 0;
            if (k > 1 && i + copies*k < n) {
                const size_t max_partial = std::min(k - 1, n - (i + copies*k));
                while (partial < max_partial &&
                       canonical(seq[i + copies*k + partial]) &&
                       std::toupper(static_cast<unsigned char>(seq[i + partial])) ==
                           std::toupper(static_cast<unsigned char>(seq[i + copies*k + partial]))) {
                    ++partial;
                }
            }
            const size_t span = copies*k + partial;
            const bool complete = copies >= MIN_NUM_REPEATS &&
                                  copies*k >= MS_MIN_NUM_BASES;
            const bool fractional = partial > 0 &&
                                    copies >= MIN_FULL_COPIES_WITH_PARTIAL &&
                                    span >= MIN_PARTIAL_REPEAT_BASES;
            if ((!complete && !fractional) ||
                copies > std::numeric_limits<uint16_t>::max()) continue;
            if (span > best_length || (span == best_length && k < best.motif.size())) {
                best.is_valid = true;
                best.start = static_cast<uint32_t>(i);
                best.end = static_cast<uint32_t>(i + span - 1);
                best.motif = std::string(seq+i, k);
                best.n_repeats = static_cast<uint16_t>(copies);
                best.terminal_partial_bases = static_cast<uint8_t>(partial);
                best_length = span;
            }
        }
        if (best.is_valid) {
            out.push_back(best);
            i += best_length;
        } else {
            ++i;
        }
    }
}

} // namespace

TandemRepeat find_microsatellite_in_reference(const char* seq, size_t n) {
    std::vector<TandemRepeat> repeats;
    scan_sequence(seq, n, repeats);
    return repeats.empty() ? TandemRepeat{} : repeats.front();
}

void find_microsatellites_in_reference(
    const char* seq, size_t seq_end, size_t seq_offset,
    std::vector<TandemRepeat>& mss) {
    std::vector<TandemRepeat> repeats;
    scan_sequence(seq, seq_end, repeats);
    for (auto& tr : repeats) {
        tr.start += static_cast<uint32_t>(seq_offset);
        tr.end += static_cast<uint32_t>(seq_offset);
        mss.push_back(std::move(tr));
    }
}

void find_microsatellites(
    const std::vector<std::pair<uint32_t, BaseInfo>>& bases,
    std::vector<TandemRepeat>& mss) {
    // Only contiguous aligned, non-insertion/non-deletion reference coordinates
    // may be interpreted as an exact genomic tandem repeat. Indels and missing
    // coverage split the run; never close the gap by concatenating bases.
    std::string run;
    uint32_t run_start = 0;
    uint32_t previous = 0;
    auto flush = [&]() {
        if (!run.empty()) {
            find_microsatellites_in_reference(run.data(), run.size(), run_start, mss);
            run.clear();
        }
    };
    for (const auto& item : bases) {
        const uint32_t pos = item.first;
        const BaseInfo& b = item.second;
        if (b.deleted || b.inserted || b.soft_clipped || !canonical(b.c)) {
            flush();
            previous = 0;
            continue;
        }
        if (!run.empty() && (pos <= previous || pos - previous != 1)) flush();
        if (run.empty()) run_start = pos;
        run.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(b.c))));
        previous = pos;
    }
    flush();
}

void find_microsatellites(
    const std::vector<std::pair<uint32_t, BaseInfo>>& b1,
    const std::vector<std::pair<uint32_t, BaseInfo>>& b2,
    std::vector<TandemRepeat>& mss) {
    // Do not concatenate mates across an unsequenced template gap.
    find_microsatellites(b1, mss);
    find_microsatellites(b2, mss);
    std::sort(mss.begin(), mss.end(), [](const TandemRepeat& a, const TandemRepeat& b) {
        if (a.start != b.start) return a.start < b.start;
        if (a.end != b.end) return a.end < b.end;
        return a.motif < b.motif;
    });
    mss.erase(std::unique(mss.begin(), mss.end(), [](const TandemRepeat& a, const TandemRepeat& b) {
        return a.start == b.start && a.end == b.end && a.motif == b.motif;
    }), mss.end());
}
