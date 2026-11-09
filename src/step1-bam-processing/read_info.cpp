//  read_info.cpp
//
//  Copyright 2024 Antoine Passemiers <antoine.passemiers@gmail.com>
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

#include "read_info.h"

#include <cstring>


void ReadInfo::trim_alignment_inplace_(const Parameters& params, bool left) {
    bool uninformative = true;
    size_t lookahead = params.trim_lookahead;
    size_t n_invalid_bases = 0;
    for (size_t i = 0; i < bases.size(); i++) {
        size_t j = (left) ? i : bases.size() - i - 1;
        BaseInfo& base = bases.at(j).second;
        const BaseInfo& ref_base = ref_bases.at(i).second;
        if (base.uninformative || base.putative_snp || (base.c == 'N') || base.inserted || base.deleted) {
            uninformative = true;
        } else if (params.c_conversion && (base.is_chg_site() || base.is_chh_site())) {
            if (base.methylated) uninformative = true;  // Methylation in CHH/CHG contexts is unlikely
            else uninformative = false;
        } else if (params.c_conversion && base.is_cpg_site()) {
            uninformative = false;
        } else if (base.c != ref_base.c) {
            uninformative = true;
        }
        if (uninformative) {
            base.uninformative = true;
            lookahead = params.trim_lookahead;
            n_invalid_bases = i + 1;
        } else {
            lookahead--;
            if (lookahead == 0) break;
        }
        if (++i > params.trim_lookahead) break;
    }

    for (size_t i = 0; i < n_invalid_bases; i++) {
        size_t j = (left) ? i : bases.size() - i - 1;
        BaseInfo &base = bases.at(j).second;
        base.uninformative = true;
    }
}


void ReadInfo::trim_alignment_inplace(const Parameters& params) {
    trim_alignment_inplace_(params, true);
    trim_alignment_inplace_(params, false);
}
