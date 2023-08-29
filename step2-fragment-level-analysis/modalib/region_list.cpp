//  region_list.cpp
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

#include "region_list.h"
#include "utils.h"

#include <string>


namespace modalib {

void _RegionList::add_region(uint8_t chr_id, uint32_t start, uint32_t end) {
    regions.push_back(Region(chr_id, start, end));
}

void _RegionList::load(const std::string &filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filepath);
    } else {
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;
            std::size_t i = 0;
            std::string chr_name;
            uint32_t start;
            uint32_t end;
            while(std::getline(ss, token, ',')) {
                switch (i) {
                    case 0:
                        chr_name = token;
                        break;
                    case 1:
                        start = std::stoi(token);
                        break;
                    case 2:
                        end = std::stoi(token);
                        break;
                }
                i++;
            }
            add_region(chr_name, start, end);
        }
    }
}

bool _RegionList::contains(const FragmentInfo &frag_info) {
    if (end_reached()) return false;
    int res;
    while ((cursor < regions.size()) && ((res = regions[cursor].cmp(frag_info)) == -1)) {
        cursor++;
    }
    if (cursor < regions.size()) {
        return (res != 1) && (res != -1);
    } else {
        return false;
    }
}

Region& _RegionList::get_current_region() {
    return regions[cursor];
}

} // namespace modalib
