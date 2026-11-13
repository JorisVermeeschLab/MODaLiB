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
#include <regex>
#include <stdexcept>

#include "htslib/sam.h"
#include "htslib/faidx.h"

#include "strand.h"


struct Parameters {

    // Minimum MAPQ value to include a read
    int min_mapq = 10;

    // Minimum Phred score to include a base
    int min_phred = 5;

    // Whether BAM file contains single-read data (instead of paired-end)
    int single_read = 0;

    // Whether to call methylation from converted cytosines.
    int c_conversion = 0;

    // Whether to call methylation from the MM and ML auxiliary fields.
    int use_mm_ml = 0;

    // If `use_mm_ml` is true, this threshold defines the minimum values in
    // the ML auxiliary field for a base to be considered as modified (i.e., methylated).
    uint8_t ml_cutoff = 168;  // 0.66 * 255

    // Maximum fragment size
    int max_fragment_size = 1000;

    // Read group filter
    bool use_read_group = false;
    std::regex read_group;

    // Input files
    char* bam_filename = nullptr;
    htsFile* fp = nullptr;
    hts_idx_t* bai = nullptr;
    char* fasta_filename = nullptr;

    // Output files
    FILE** output_fp = nullptr;
    FILE** c_output_fp = nullptr;
    bool write_cytosines = false;
    FILE* msi_output_fp = nullptr;
    bool write_msi = false;
    uint32_t msi_flank = 3;

    // Inclusion bounds for the 4 types of (converted) strands.
    // Useful resources:
    // 1) https://software.broadinstitute.org/software/igv/interpreting_bisulfite_mode
    // 2) https://www.strand-ngs.com/files/manual/reference/appx_methyl_seq.html#SECTION00530000000000000000
    int ot_bounds[4] = {0};
    int ob_bounds[4] = {0};
    int ctot_bounds[4] = {0};
    int ctob_bounds[4] = {0};

    // Parameters for the algorithmic trimming of the reads
    size_t trim_lookahead = 2;
    size_t trim_max_n_bases = 12;

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
    void parse_ctot_bounds(const char* arg) { parse_bounds(arg, ctot_bounds); }
    void parse_ctob_bounds(const char* arg) { parse_bounds(arg, ctob_bounds); }
    static void parse_bounds(const char* arg, int* bounds) {
        // Reject more than four values *before* writing beyond bounds[4].
        size_t count = 0;
        std::string remaining(arg);
        size_t pos = 0;
        while (true) {
            pos = remaining.find(',');
            const std::string token = remaining.substr(0, pos);
            if (token.empty() || count == 4) {
                throw std::invalid_argument("Expected exactly four nonempty inclusion bounds");
            }
            bounds[count++] = std::stoi(token);
            if (pos == std::string::npos) break;
            remaining.erase(0, pos + 1);
        }
        if (count != 4) {
            throw std::invalid_argument("Expected exactly four inclusion bounds");
        }
    }

};

#endif // PARAMS_H_
