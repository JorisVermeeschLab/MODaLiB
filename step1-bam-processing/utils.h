//  methylation.h
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

#ifndef UTILS_H_
#define UTILS_H_

#include "htslib/sam.h"
#include "htslib/faidx.h"

#include <sstream>
#include <string>
#include <utility>


/*! @function
 @abstract Determine the strand from which a read arose
 @param	b	Pointer to an alignment
 @returns	1, original top; 2, original bottom; 3, complementary to the original top; 4, complementary to the original bottom
 @discussion There are two methods used to determine the strand of origin. The 
 first method uses XG auxiliary tags as output by Bison and Bismark. For details
 on how strand is determined from this, see the source code for this function or
 the documentation for either Bison or bismark. This method can handle non-
 directional libraries. If an XG auxiliary tag is not present, then the
 orientation of each alignment (and whethers it's read #1 or #2 of a pair, if
 applicable) determines the strand.
*/
int get_strand(bam1_t *b) {
    char *XG = (char *) bam_aux_get(b, "XG");
    //Only bismark uses the XG tag like this. Some other aligners use it for other purposes...
    if(XG != NULL && *(XG+1) != 'C' && *(XG+1) != 'G') XG = NULL;
    if(XG == NULL) { //Can't handle non-directional libraries!
        if(b->core.flag & BAM_FPAIRED) {
            if((b->core.flag & 0x50) == 0x50) return 2; //Read1, reverse comp. == OB
            else if(b->core.flag & 0x40) return 1; //Read1, forward == OT
            else if((b->core.flag & 0x90) == 0x90) return 1; //Read2, reverse comp. == OT
            else if(b->core.flag & 0x80) return 2; //Read2, forward == OB
            return 0; //One of the above should be set!
        } else {
            if(b->core.flag & 0x10) return 2; //Reverse comp. == OB
            return 1; //OT
        }
    } else {
        if(*(XG+1) == 'C') { //OT or CTOT, due to C->T converted genome
            if((b->core.flag & 0x51) == 0x41) return 1; //Read#1 forward == OT
            else if((b->core.flag & 0x51) == 0x51) return 3; //Read #1 reverse == CTOT
            else if((b->core.flag & 0x91) == 0x81) return 3; //Read #2 forward == CTOT
            else if((b->core.flag & 0x91) == 0x91) return 1; //Read #2 reverse == OT
            else if(b->core.flag & 0x10) return 3; //Single-end reverse == CTOT
            else return 1; //Single-end forward == OT
        } else {
            if((b->core.flag & 0x51) == 0x41) return 4; //Read#1 forward == CTOB
            else if((b->core.flag & 0x51) == 0x51) return 2; //Read #1 reverse == OB
            else if((b->core.flag & 0x91) == 0x81) return 2; //Read #2 forward == OB
            else if((b->core.flag & 0x91) == 0x91) return 4; //Read #2 reverse == CTOB
            else if(b->core.flag & 0x10) return 2; //Single-end reverse == OB
            else return 4; //Single-end forward == CTOB
        }
    }
}


std::string get_read_as_string(bam1_t *b) {
	uint8_t* seq = bam_get_seq(b);
	int seq_len = b->core.l_qseq;

	std::stringstream ss;
	for (int i = 0; i < seq_len; i++) {
		int base = bam_seqi(seq, i);
		switch (base) {
		case 2:
			ss << "C";
			break;
		case 8:
			ss << "T";
			break;
		case 4:
			ss << "G";
			break;
		case 1:
			ss << "A";
			break;
		default:
			ss << "N";
		}
	}
	return ss.str();
}


#endif // UTILS_H_