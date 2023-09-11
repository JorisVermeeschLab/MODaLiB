#ifndef UTILS_H_
#define UTILS_H_

#include "htslib/sam.h"
#include "htslib/faidx.h"

#include <sstream>
#include <string>
#include <utility>


const int STRAND_OT = 1;
const int STRAND_OB = 2;
const int STRAND_CTOT = 3;
const int STRAND_CTOB = 4;


/*! @function
 *
 * If the library is directional, only reads originating from the original top strand (OT)
 * or the original bottom strand (OB) (up to a few cytosine conversions) will be sequenced.
 *
 @abstract Determine the strand from which a read arose
 @param	b	Pointer to an alignment
 @returns	1, original top; 2, original bottom; 3, complementary to the original top; 4, complementary to the original bottom
 @discussion There are two methods used to determine the strand of origin. The 
 first method uses XG auxiliary tags as output by Bison and Bismark. For details
 on how strand is determined from this, see the source code for this function or
 the documentation for either Bison or bismark. This method can handle non-
 directional libraries. If an XG auxiliary tag is not present, then the
 orientation of each alignment (and whether it's read #1 or #2 of a pair, if
 applicable) determines the strand.
*/
int get_strand(bam1_t *b) {
    uint16_t flag = b->core.flag;
    char* XG = (char *) bam_aux_get(b, "XG");
    //Only bismark uses the XG tag like this. Some other aligners use it for other purposes...
    if (XG != nullptr && *(XG+1) != 'C' && *(XG+1) != 'G') XG = nullptr;
    if (XG == nullptr) {
        // This code block does not work for non-directional libraries
        if (flag & BAM_FPAIRED) {  // paired-end
            if (flag & BAM_FREAD1) {  // Read mapped to the reverse strand
                if (flag & BAM_FREVERSE) return STRAND_OB;
                else return STRAND_OT;
            } else if (flag & BAM_FREAD2) {
                if (flag & BAM_FREVERSE) return STRAND_OT;
                else return STRAND_OB;
            } else {
                return STRAND_OT;  // TODO: do something smarter
            }
        } else {  // Single-end
            if (flag & BAM_FREVERSE) return STRAND_OB;
            return STRAND_OT;
        }
    } else {
        if (*(XG+1) == 'C') {
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
        } else {
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


std::string to_string(bam1_t *b) {
	uint8_t* seq = bam_get_seq(b);
	int seq_len = b->core.l_qseq;

	std::stringstream ss;
	for (int i = 0; i < seq_len; i++) {
		ss << seq_nt16_str[bam_seqi(seq, i)];
	}
	return ss.str();
}


#endif // UTILS_H_