//  incremental.h
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

#ifndef MODALIB_INCREMENTAL_H
#define MODALIB_INCREMENTAL_H

#include <deque>

#include "utils.h"
#include "fragment_info.h"

namespace modalib {

class IncrementalMean {
private:
    float running_sum = 0.0f;
    float total_weight = 0.0f;
    size_t max_window_size;
    std::deque<std::pair<FragmentInfo, std::pair<float, float>>> window;
public:
    explicit IncrementalMean(size_t max_window_size) : running_sum(0.0f), total_weight(0.0f), max_window_size(max_window_size), window() {}
    void add(const FragmentInfo& frag, float value, float weight) {
        window.push_back({frag, {value, weight}});
        running_sum += value;
        total_weight += weight;
        update();
    }
    bool contains(const FragmentInfo& frag) const {
        if (window.empty()) return true;
        const FragmentInfo& newest = window.back().first;
        if (newest.start_chr_id != frag.start_chr_id) return false;
        long diff = static_cast<long>(newest.start_pos) - static_cast<long>(frag.start_pos);
        if (std::labs(diff) >= max_window_size) return false;
        return true;
    }
    void update() {
        while (!contains(window.front().first)) {
            auto& pair = window.front();
            running_sum -= pair.second.first;
            total_weight -= pair.second.second;
            window.pop_front();
        }
        if (window.empty()) {
            running_sum = 0.0f;
            total_weight = 0.0f;
        }
    }
    inline float as_float() const {
        size_t n = window.size();
        if (n == 0) return 0.0f;
        return running_sum / static_cast<float>(n);
    }
    inline float get_weight() const { return total_weight; }
    inline size_t size() const { return window.size(); }
};

} // namespace modalib

#endif //MODALIB_INCREMENTAL_H
