//  strand.h
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

#ifndef MODALIB_STRAND_H
#define MODALIB_STRAND_H

#include <cstring>
#include <sstream>
#include <string>
#include <utility>

#include "htslib/sam.h"
#include "htslib/faidx.h"


const int STRAND_OT = 1;
const int STRAND_OB = 2;
const int STRAND_CTOT = 3;
const int STRAND_CTOB = 4;


int get_bismark_conversion_from_tag(const char* tag);


/**
 * Determine the strand of origin of a read,
 * assuming Bismark was used to align the reads.
 *
 * XR == "CT" denotes that the read has been aligned in its original form.
 * XR == "GA" denotes that the reverse complement has been aligned instead.
 * XG == "CT" denotes that the read has been aligned to the OT of the reference.
 * XG == "GA" denotes that the read has been aligned to the OB of the reference.
 * Using these 2 tags jointly is sufficient to disentable between OT, OB, CTOT and CTOB.
 *
 * @param b Alignment for current read.
 * @return The strand of origin, if optional fields from
 *      Bismark have been found, -1 otherwise.
 */
int get_strand_using_bismark(const bam1_t* b);


/**
 * Determine the strand of origin of a read.
 *
 * `BAM_FREVERSE` indicates whether SEQ has been reverse complemented and QUAL reversed.
 * As a reminder, in Illumina sequencing the reverse strand (`BAM_FREVERSE`, not to be confused
 * with the bottom strand -- OB) is obtained by bridge amplification of the forward strand
 * (not to be confused with the top strand -- OT).
 * Therefore, if the read is reversed (`BAM_FREVERSE`) and coming from the forward
 * strand (`BAM_FREAD1`), then it means that the read actually originates from the
 * bottom strand (OB).
 *
 * Schematically, we have:
 * OT/OB/CTOT/CTOB --[bridge amplification]--> R1/R2 --[alignment]--> forward/reverse
 *
 * @param b Alignment for current read.
 * @return The strand of origin.
 */
int get_strand(const bam1_t *b);


std::string to_string(bam1_t *b);


#endif // MODALIB_STRAND_H
