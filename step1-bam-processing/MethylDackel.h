// IMPORTANT NOTE: this file is a modified version of the header file from the MethylDackel software

#include <inttypes.h>
#include <stdio.h>
#include <zlib.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"

#include "utils.h"

#if HTS_VERSION < 101100
#error "The minimum supported version of htslib is 1.11!"
#endif


//These are needed to handle multiple threads
extern pthread_mutex_t positionMutex;
extern pthread_mutex_t outputMutex;
extern uint32_t globalTid;
extern uint32_t globalPos;
extern uint32_t globalEnd;
extern uint32_t bin;
extern uint32_t outputBin;

/*! @typedef
 @abstract Global configuration information structure
 @field	minMapq	Minimum MAPQ value to include a read (-q)
 @field	minPhred	Minimum Phred score to include a base (-p)
 @field	minDepth	Minimum depth for outputing methylation metrics
 @field ignoreFlags     Mask that's logically &ed with and ignored if > 0. Defaults to 0xF00.
 @field requireFlags    Mask that's logically &ed with and ignored if < mask. Defaults to 0, which means ignore.
 @field output_fp	Output file pointers (to CpG, CHG, and CHH, respectively)
 @field	reg	A region specified by -r
 @field BAMName	The BAM file name
 @field fp	Input file pointer (must be a BAM or CRAM file)
 @field	bai	The index for fp
 @field FastaName	The fasta file name
 @field bounds	Trimming bounds
 @field absoluteBounds	Absolute trimming bounds
 @field nThreads	Number of threads in use.
 @field chunkSize	The number of bases processed by each thread at a time (can be adjusted a bit to ensure CpGs/CHGs aren't split between processors)
*/
typedef struct {
    int minMapq = 10;
    int minPhred = 5;
    int minDepth;
    int ignoreFlags = 0;
    int requireFlags = 0;
    int singleRead = 0;
    int use_mm_ml = 0;
    uint8_t ml_cutoff = 1;
    FILE **output_fp;
    char *reg = NULL;
    char *BAMName;
    htsFile *fp = NULL;
    hts_idx_t *bai = NULL;
    char *FastaName;
    int bounds[16] = {0};
    int absoluteBounds[16] = {0};
    int nThreads = 1;
    unsigned long chunkSize = 1000000;
} Config;

//Used internally to parse things like --OT 0,1,2,3
void parseBounds(char *s2, int *vals, int mult);

//returns <0 if region0 comes before region1
//returns >0 if region0 comes after region1
//returns 0 on overlap
static inline int64_t compareRegions(int32_t tid0, int32_t start0, int32_t end0, int32_t tid1, int32_t start1, int32_t end1) {
    if(tid0 != tid1) return ((int64_t) tid0)-((int64_t) tid1);
    if(start0 < start1 && end0 >= start1) return 0;
    if(start0 >= start1 && start0 < end1) return 0;
    return ((int64_t) start0)-((int64_t) start1);
}

bam1_t *trimAlignmentSimple(bam1_t *b, int bounds[16]) {
    int strand = get_strand(b)-1;
    int i, lb, rb;
    uint8_t *qual = bam_get_qual(b);

    if(b->core.flag & BAM_FREAD2) {
        lb = bounds[4*strand+2];
        rb = bounds[4*strand+3];
    } else {
        lb = bounds[4*strand];
        rb = bounds[4*strand+1];
    }

    lb = (lb<b->core.l_qseq) ? lb : b->core.l_qseq;

    if (lb) for(i=0; i<lb; i++) qual[i] = 0;
    if (rb) for(i=rb; i<b->core.l_qseq; i++) qual[i] = 0;
    return b;
}

void parseBounds(char *s2, int *vals, int mult) {
    char *p, *s = strdup(s2), *end;
    int i, v;
    long tempV;

    p = strtok(s, ",");
    tempV = strtol(p, &end, 10);
    if((errno == ERANGE && (tempV == LONG_MAX || tempV == LONG_MIN)) || (errno != 0 && tempV == 0) || end == p) v = -1;
    else if(tempV > INT_MAX || tempV < LONG_MIN) v = -1;
    else v = tempV;

    if(v>=0) vals[4*mult] = v;
    else if (v == -1) vals[4*mult] = 1;
    else {
        fprintf(stderr, "Invalid bounds string, %s\n", s2);
        free(s);
        return;
    }
    for(i=1; i<4; i++) {
        p = strtok(NULL, ",");
        tempV = strtol(p, &end, 10);
        if((errno == ERANGE && (tempV == LONG_MAX || tempV == LONG_MIN)) || (errno != 0 && tempV == 0) || end == p) v = -1;
        else if(tempV > INT_MAX || tempV < LONG_MIN) v = -1;
        else v = tempV;

        if(v>=0) vals[4*mult+i] = v;
        else if (v == -1) {
            if (i == 2) vals[4*mult+i] = 0;
            else vals[4*mult+i] = 0;
        }
        else {
            fprintf(stderr, "Invalid bounds string, %s\n", s2);
            free(s);
            return;
        }
    }
    free(s);
}
