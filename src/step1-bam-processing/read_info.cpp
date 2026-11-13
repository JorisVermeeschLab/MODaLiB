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
    // Base/ref-base vectors are aligned by index by process_read().
    const size_t n = std::min(bases.size(), ref_bases.size());
    if (!n || !params.trim_lookahead || !params.trim_max_n_bases) return;

    const size_t limit = std::min(n, params.trim_max_n_bases);
    size_t last_invalid = 0;
    size_t consecutive_valid = 0;
    for (size_t i = 0; i < limit; ++i) {
        const size_t j = left ? i : n - 1 - i;
        const BaseInfo& base = bases[j].second;
        const BaseInfo& ref_base = ref_bases[j].second;
        bool invalid = base.uninformative || base.putative_snp ||
                       base.c == 'N' || base.inserted || base.deleted;
        if (!invalid && params.c_conversion &&
            (base.is_chg_site() || base.is_chh_site())) {
            invalid = base.methylated;
        } else if (!invalid && params.c_conversion && base.is_cpg_site()) {
            invalid = false;
        } else if (!invalid) {
            invalid = base.c != ref_base.c;
        }
        if (invalid) {
            last_invalid = i + 1;
            consecutive_valid = 0;
        } else if (++consecutive_valid >= params.trim_lookahead) {
            break;
        }
    }
    for (size_t i = 0; i < last_invalid; ++i) {
        const size_t j = left ? i : n - 1 - i;
        bases[j].second.uninformative = true;
    }
}

void ReadInfo::trim_alignment_inplace(const Parameters& params) {
    trim_alignment_inplace_(params, true);
    trim_alignment_inplace_(params, false);
}
