//  counter.h
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

#ifndef MODALIB_COUNTER_H
#define MODALIB_COUNTER_H

#include <algorithm>
#include <utility>
#include <unordered_map>
#include <deque>
#include <vector>
#include <iostream>

#include "fragment_info.h"
#include "bounds.h"

namespace modalib {

const long WINDOW_SIZE = 50000;

class Counter {
private:
    long total = 0;
    std::unordered_map<int, std::pair<long, long>> bounds;

    std::vector<long> counts;
    std::vector<long> total_counts;
    std::deque<FragmentInfo> window;
    static inline int gc_index(const FragmentInfo& frag) {
        const int gc = static_cast<int>(std::roundf(100.0f * frag.get_gc_ratio()));
        return std::max(0, std::min(gc, 100));
    }
    bool contains(const FragmentInfo& frag) {
        if (window.size() == 0) return true;
        const FragmentInfo& newest = window.back();
        if (newest.start_chr_id != frag.start_chr_id) return false;
        if (std::labs(newest.start_pos - frag.start_pos) >= WINDOW_SIZE) return false;
        return true;
    }
    void update() {
        while ((window.size() > 0) && (!contains(window.front()))) {
            const int gc = gc_index(window.front());
            counts[gc]--;
            window.pop_front();
        }
    }
public:
    Counter() : total(0), bounds(), counts(101, 0), total_counts(101, 0), window() {}
    Counter(const Counter& other) : total(other.total), bounds(other.bounds), counts(other.counts), total_counts(other.total_counts), window(other.window) {}
    void add(const FragmentInfo& frag) {
        const int gc = gc_index(frag);
        counts[gc]++;
        window.push_back(frag);
    } // TODO: difference between add and update?
    void inc_total_counts(const FragmentInfo& frag) {
        const int gc = gc_index(frag);
        total_counts[gc]++;

        const int partition_id = coords_to_partition(frag.start_chr_id, frag.start_pos, NUM_PARTITIONS);
        auto it = bounds.find(partition_id);
        if (it != bounds.end()) {
            it->second.second = total + 1;
        } else {
            bounds[partition_id] = std::pair<long, long>(total, total + 1);
        }

        total++;
    }
    float get_density(const FragmentInfo& frag) noexcept {
        const int gc = gc_index(frag);
        if (total_counts[gc] == 0) return 0.0;
        update();
        // std::cout << counts[gc] << ", " << total_counts[gc] << ", ";
        return static_cast<float>(counts[gc]) / static_cast<float>(total_counts[gc]);
    }
    std::pair<long, long> get_bounds(int partition_id) {
        auto it = bounds.find(partition_id);
        if (it != bounds.end()) {
            return it->second;
        } else {
            return std::pair<long, long>(0, 0);
        }
    }
};

} // namespace modalib

#endif  // MODALIB_COUNTER_H
