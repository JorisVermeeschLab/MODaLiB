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

struct BaseInfo {

    // Char representation of the base.
    // If cytosine conversion has occured (in the context of methylation sequencing),
    // then the base has been reverted back to the correct base (in order to match
    // the reference genome).
    char c = 'N';

    // Phred score
    uint8_t phred_score = 0;

    // Whether base corresponds to the C of a CpG site
    uint8_t is_cpg_site : 1;

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
        is_cpg_site = 0;
        methylated = 0;
        inserted = 0;
        deleted = 0;
        putative_snp = 0;
        uninformative = 0;
        soft_clipped = 0;
    }

    inline bool is_c_or_g() const { return (c == 'C') or (c == 'G'); }
    inline char type_as_char() const {
        if (soft_clipped) {
            return 'S';
        } else if (inserted) {
            return '/';
        } else if (deleted) {
            return '\\';
        } else if (is_cpg_site) {
            if (methylated) return 'M';
            else return 'U';
        } else if (putative_snp) {
            return 'X';
        } else if (uninformative) {
            return '0';
        } else {
            return ' ';
        }
    }
};

#endif //MODALIB_BASE_H
