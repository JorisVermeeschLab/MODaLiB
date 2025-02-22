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

#include <algorithm>
#include <string>
#include <stdexcept>


void MSS::add_region(uint8_t chr_id, uint32_t start, uint32_t end, const std::string& motif, size_t n_repeats) {
    int new_region_id = regions.size();
    regions.emplace_back(new_region_id, chr_id, start, end, motif, n_repeats);
    if (region_ids_by_chr.find(chr_id) == region_ids_by_chr.end()) {
        region_ids_by_chr[chr_id] = { new_region_id };
    } else {
        region_ids_by_chr[chr_id].push_back(new_region_id);
    }
    cursors[chr_id] = 0;
}

void MSS::load(const std::string &filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filepath);
    } else {
        uint8_t last_chr_id = 0;
        uint32_t last_start = 0;
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;
            std::size_t i = 0;
            uint8_t chr_id = 0;
            uint32_t start = 0;
            uint32_t end = 0;
            std::string motif;
            while (std::getline(ss, token, '\t')) {
                switch (i) {
                    case 0:
                        chr_id = parse_chr_name(token.c_str());
                        break;
                    case 1:
                        start = std::stoi(token);
                        break;
                    case 2:
                        end = std::stoi(token);
                        break;
                    case 3:
                        motif = token;
                        break;
                }
                i++;
            }

            auto p = parse_mss_motif(motif);
            add_region(chr_id, start, end, p.first, p.second);

            // Check if regions are properly ordered
            bool sorted = true;
            if (chr_id < last_chr_id) {
                sorted = false;
            } else if (chr_id > last_chr_id) {
                sorted = true;
            } else {
                sorted = (last_start <= start);
            }
            if (!sorted) {
                throw std::runtime_error("Region list is not sorted");
            }
        }
    }
    std::cout << "Loaded reference MS list." << std::endl;
}

void MSS::reset() {
    for (auto& it: cursors) it.second = 0;
}

Overlap MSS::overlap(const Fragment& frag_info, int region_id) {
    if (region_id >= regions.size()) return Overlap();
    return regions[region_id].compute_overlap(frag_info);
}

Overlap MSS::search_overlap(const Fragment &frag_info) {
    Overlap res;
    res.found_contig = true;
    res.reached_contig_end = false;

    // Ensure contig is in the list of regions
    if (region_ids_by_chr.find(frag_info.chr_id) == region_ids_by_chr.end()) {
        res.found_contig = false;
        return res;
    }
    if (cursors.find(frag_info.chr_id) == cursors.end()) {
        res.found_contig = false;
        return res;
    }
    std::vector<int>& region_ids = region_ids_by_chr[frag_info.chr_id];

    // Ensure there is at least one region in that contig
    if (region_ids.empty()) {
        res.reached_contig_end = true;
        return res;
    }

    // Ensure cursor is not out-of-bounds
    int& cursor = cursors[frag_info.chr_id];
    cursor = std::min(cursor, (int) region_ids.size() - 1);

    // Ensure that left end is not left (outside) to current region
    while ((cursor > 0) && ((res = regions[region_ids[cursor]].compute_overlap(frag_info)).may_overlap_previous_region)) {
        cursor--;
    }

    // Search for the leftmost region containing the left end of the fragment
    while ((cursor < region_ids.size()) && ((res = regions[region_ids[cursor]].compute_overlap(frag_info)).rpos_type == RPOS_RIGHT_NO_OVERLAP)) {
        cursor++;
    }

    if (cursor >= region_ids.size()) {
        res = Overlap();
        res.found_contig = true;
        res.reached_contig_end = true;
    }
    return res;
}
