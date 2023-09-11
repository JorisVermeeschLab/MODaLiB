//  params.h
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

#ifndef PARAMS_H_
#define PARAMS_H_

#include <algorithm>
#include <cinttypes>
#include <cstdio>
#include <zlib.h>
#include <iostream>
#include <string>
#include <stdexcept>

#include "htslib/sam.h"
#include "htslib/faidx.h"

#include "utils.h"


struct Parameters {

    // Minimum MAPQ value to include a read
    int min_mapq = 10;

    // Minimum Phred score to include a base
    int min_phred = 5;

    // Whether BAM file contains single-read data (instead of paired-end)
    int single_read = 0;

    // Whether to call methylation from the MM and ML auxiliary fields.
    // Otherwise, methylation is called from the converted cytosines.
    int use_mm_ml = 0;

    // If `use_mm_ml` is true, this threshold defines the minimum values in
    // the ML auxiliary field for a base to be considered as modified (i.e., methylated).
    uint8_t ml_cutoff = 168;  // 0.66 * 255

    FILE** output_fp = nullptr;
    char* bam_filename = nullptr;
    htsFile* fp = nullptr;
    hts_idx_t* bai = nullptr;
    char* fasta_filename = nullptr;

    // Inclusion bounds for the 4 types of (converted) strands.
    // The difference between the 4 types of strands are expained is illustrated here:
    // https://software.broadinstitute.org/software/igv/interpreting_bisulfite_mode
    int ot_bounds[4] = {0};
    int ob_bounds[4] = {0};
    int ctot_bounds[4] = {0};
    int ctob_bounds[4] = {0};

    int n_threads = 1;
    unsigned long partition_size = 1000000;

    const int* get_bounds(int strand) const {
        const int* bounds;
        switch (strand) {
            case STRAND_OT:  // Original top strand
                bounds = ot_bounds;
                break;
            case STRAND_OB:  // Original bottom strand
                bounds = ob_bounds;
                break;
            case STRAND_CTOT:  // Complementary to the top strand
                bounds = ctot_bounds;
                break;
            case STRAND_CTOB:  // Complementary to the bottom strand
                bounds = ctob_bounds;
                break;
            default:
                bounds = ot_bounds;  // TODO: throw error
        }
        return bounds;
    }
    void parse_ot_bounds(const char* arg) { parse_bounds(arg, ot_bounds); }
    void parse_ob_bounds(const char* arg) { parse_bounds(arg, ob_bounds); }
    void parse_ctot_bounds(const char* arg) { parse_bounds(arg, ctob_bounds); }
    void parse_ctob_bounds(const char* arg) { parse_bounds(arg, ctot_bounds); }
    static void parse_bounds(const char* arg, int* bounds) {
        size_t i = 0;
        std::string s(arg);
        size_t pos;
        std::string token;
        while ((pos = s.find(',')) != std::string::npos) {
            token = s.substr(0, pos);
            bounds[i] = std::stoi(token);
            s.erase(0, pos + 1);
            i++;
        }
        if (s.length() > 0) {
            bounds[i++] = std::stoi(s);
        }
        if (i != 4) {
            throw std::invalid_argument("Inclusion bounds parameter expected 4 integers. Got " + std::to_string(i) + ".");
        }
    }
};

#endif // PARAMS_H_
