//  region_list.h
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

#ifndef MODALIB_REGION_LIST_H
#define MODALIB_REGION_LIST_H

#include "fragment_info.h"
#include "site.h"
#include "overlap.h"

#include <fstream>
#include <unordered_map>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

namespace modalib {

class _RegionList {
private:
    std::vector<Region> regions;
    std::unordered_map<uint8_t, std::vector<int>> region_ids_by_chr;
    std::unordered_map<uint8_t, int> cursors;
public:
    _RegionList(): regions(), region_ids_by_chr(), cursors() {}
    _RegionList(const _RegionList& other) = default;
    void load(const std::string &filepath);
    void add_region(uint8_t chr_id, uint32_t start, uint32_t end);
    Overlap search_overlap(const FragmentInfo&);
    Overlap overlap(const FragmentInfo&, int region_id);
    void reset();
    inline const Region& at(size_t region_id) const { return regions.at(region_id); }
    inline std::size_t get_num_regions() const { return regions.size(); }
};

} // namespace modalib

#endif // MODALIB_REGION_LIST_H
