//  strand.cpp
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

#include "strand.h"

#include <cstring>
#include <sstream>
#include <string>
#include <utility>

#include "htslib/sam.h"
#include "htslib/faidx.h"


int get_bismark_conversion_from_tag(const char* tag) {
    if (tag == nullptr) return 0;
    if (std::strlen(tag) >= 3) {
        if (((*(tag + 1)) == 'C') && ((*(tag + 2)) == 'T')) {
            return 1;
        } else if (((*(tag + 1)) == 'G') && ((*(tag + 2)) == 'A')) {
            return -1;
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}


int get_strand_using_bismark(const bam1_t* b) {
    uint16_t flag = b->core.flag;

    const char* XG = (char *) bam_aux_get(b, "XG");
    int genome_conversion_state = get_bismark_conversion_from_tag(XG);
    if (genome_conversion_state == 0) return -1;
    const char* XR = (char *) bam_aux_get(b, "XR");
    int read_conversion_state = get_bismark_conversion_from_tag(XR);

    if (genome_conversion_state > 0) {  // XG == "CT"
        if (read_conversion_state > 0) {  // XR == "CT"
            return STRAND_OT;
        } else if (read_conversion_state < 0) {  // XR == "GA"
            return STRAND_CTOT;
        } else {  // Ambiguous case
            if (flag & BAM_FPAIRED) {  // Paired-end
                if (flag & BAM_FREAD1) {
                    if (flag & BAM_FREVERSE) return STRAND_CTOT;
                    else return STRAND_OT;
                } else if (flag & BAM_FREAD2) {
                    if (flag & BAM_FREVERSE) return STRAND_OT;
                    else return STRAND_CTOT;
                }
            } else {  // Single-end
                if (flag & BAM_FREVERSE) return STRAND_CTOT;
                else return STRAND_OT;
            }
        }
    } else {  // XG == "GA"
        if (read_conversion_state > 0) {  // XR == "CT"
            return STRAND_OB;
        } else if (read_conversion_state < 0) {  // XR == "GA"
            return STRAND_CTOB;
        } else {  // Ambiguous case
            if (flag & BAM_FPAIRED) {  // Paired-end
                if (flag & BAM_FREAD1) {
                    if (flag & BAM_FREVERSE) return STRAND_OB;
                    else return STRAND_CTOB;
                } else if (flag & BAM_FREAD2) {
                    if (flag & BAM_FREVERSE) return STRAND_CTOB;
                    else return STRAND_OB;
                }
            } else {  // Single-end
                if (flag & BAM_FREVERSE) return STRAND_OB;
                else return STRAND_CTOB;
            }
        }
    }
}


int get_strand(const bam1_t *b) {
    uint16_t flag = b->core.flag;

    // Check whether Bismark has been used to align the reads
    // and use aligner-specific optional fields (in case) to call methylation
    int strand = get_strand_using_bismark(b);
    if (strand >= 0) return strand;

    // This code block does not work for non-directional libraries
    if (flag & BAM_FPAIRED) {  // paired-end
        if (flag & BAM_FREAD1) {
            if (flag & BAM_FREVERSE) return STRAND_OB;  // Read mapped to the reverse strand
            else return STRAND_OT;
        } else if (flag & BAM_FREAD2) {
            if (flag & BAM_FREVERSE) return STRAND_OT;  // Read mapped to the reverse strand
            else return STRAND_OB;
        } else {
            return STRAND_OT;  // TODO: do something smarter
        }
    } else {  // Single-end
        if (flag & BAM_FREVERSE) return STRAND_OB;
        return STRAND_OT;
    }
}


std::string to_string(bam1_t *b) {
    uint8_t* seq = bam_get_seq(b);
    int seq_len = b->core.l_qseq;

    std::stringstream ss;
    for (int i = 0; i < seq_len; i++) {
        ss << seq_nt16_str[bam_seqi(seq, i)];
    }
    return ss.str();
}
