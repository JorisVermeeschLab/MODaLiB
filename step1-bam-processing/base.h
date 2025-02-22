//  base.h
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

#ifndef MODALIB_BASE_H
#define MODALIB_BASE_H

#include <cstdint>
#include "methylation/context.h"

struct BaseInfo {

    // Char representation of the base.
    // If cytosine conversion has occurred (in the context of methylation sequencing),
    // then the base has been reverted back to the correct base (in order to match
    // the reference genome).
    char c = 'N';

    // Phred score
    uint8_t phred_score = 0;

    // Sequence context
    Context context = NOT_C_CONTEXT;

    // Whether base is methylated
    uint8_t methylated : 1;

    // Whether base is a match (equal or different)
    uint8_t is_match : 1;

    // Whether base has been inserted
    uint8_t inserted : 1;

    // Whether base has been deleted
    uint8_t deleted : 1;

    // Whether base is a putative SNP
    uint8_t putative_snp : 1;

    // Whether base is uninformative (e.g., base is unknown, or Phred score is too low).
    // For example, whether it cannot be used to estimate the GC-content
    // of the read/fragment.
    uint8_t uninformative : 1;

    // Whether base is soft-clipped
    uint8_t soft_clipped : 1;

    BaseInfo() {
        is_match = 0;
        context = NOT_C_CONTEXT;
        methylated = 0;
        inserted = 0;
        deleted = 0;
        putative_snp = 0;
        uninformative = 0;
        soft_clipped = 0;
    }

    inline bool is_cpg_site() const { return context == CpG_CONTEXT; }
    inline bool is_chg_site() const { return context == CHG_CONTEXT; }
    inline bool is_chh_site() const { return context == CHH_CONTEXT; }
    inline bool has_cytosine_context() const { return is_cpg_site() || is_chg_site() || is_chh_site(); }
    inline bool is_c_or_g() const { return (c == 'C') or (c == 'G'); }
    inline char type_as_char() const {
        if (soft_clipped) {
            return 'S';
        } else if (inserted) {
            return '/';
        } else if (deleted) {
            return '\\';
        } else if (uninformative) {
            return '0';
        } else if (context == CpG_CONTEXT) {
            if (methylated) return 'M';
            else return 'U';
        } else if ((context == CHG_CONTEXT) || (context == CHH_CONTEXT)) {
            if (methylated) return 'm';
            else return 'u';
        } else if (putative_snp) {
            return 'X';
        } else {
            return ' ';
        }
    }
};

#endif //MODALIB_BASE_H
