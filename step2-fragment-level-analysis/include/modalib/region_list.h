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

#ifndef REGION_LIST_H_
#define REGION_LIST_H_

#include "fragment_info.h"
#include "site.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

namespace modalib {

class _RegionList {
private:
    std::vector<Region> regions;
    std::size_t cursor;
public:
    _RegionList(): regions(), cursor(0) {}
    void load(const std::string &filepath);
    void add_region(uint8_t chr_id, uint32_t start, uint32_t end);
    bool contains(const FragmentInfo&);
    std::size_t get_region_id() { return cursor; }
    Region& get_region(std::size_t i) { return regions[i]; }
    Region& get_current_region();
    void reset() { cursor = 0; }
    std::size_t get_num_regions() { return regions.size(); }
    bool end_reached() { return (cursor >= regions.size()); }
};

} // namespace modalib

#endif // REGION_LIST_H_
