#ifndef _COUNTER_H__
#define _COUNTER_H__

#include <algorithm>
#include <utility>
#include <unordered_map>
#include <deque>
#include <vector>
#include <iostream>

#include "read_info.h"
#include "bounds.h"

namespace modalib {

const long WINDOW_SIZE = 50000;

class Counter {
private:
    long total = 0;
    std::unordered_map<int, std::pair<long, long>> bounds;

    std::vector<long> counts;
    std::vector<long> total_counts;
    std::deque<ReadInfo> window;
    static inline int gc_index(const ReadInfo& read) {
        const int gc = static_cast<int>(std::roundf(100.0f * read.get_gc_ratio()));
        return std::max(0, std::min(gc, 100));
    }
    bool contains(const ReadInfo& read) {
        if (window.size() == 0) return true;
        const ReadInfo& newest = window.back();
        if (newest.start_chr_id != read.start_chr_id) return false;
        if (std::labs(newest.start_pos - read.start_pos) >= WINDOW_SIZE) return false;
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
    void add(const ReadInfo& read) {
        const int gc = gc_index(read);
        counts[gc]++;
        window.push_back(read);
    } // TODO: difference between add and update?
    void inc_total_counts(const ReadInfo& read) {
        const int gc = gc_index(read);
        total_counts[gc]++;

        const int partition_id = coords_to_partition(read.start_chr_id, read.start_pos, NUM_PARTITIONS);
        auto it = bounds.find(partition_id);
        if (it != bounds.end()) {
            it->second.second = total;
        } else {
            bounds[partition_id] = std::pair<long, long>(total, total);
        }
        total++;
    }
    float get_density(const ReadInfo& read) noexcept {
        const int gc = gc_index(read);
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

#endif // _DESERIALIZATION_H__
