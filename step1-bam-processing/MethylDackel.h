// IMPORTANT NOTE: this file is a modified version of the header file from the MethylDackel software

#include <inttypes.h>
#include <stdio.h>
#include <zlib.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"

#if HTS_VERSION < 101100
#error "The minimum supported version of htslib is 1.11!"
#endif


//These are needed to handle multiple threads
extern pthread_mutex_t positionMutex;
extern pthread_mutex_t bwMutex;
extern pthread_mutex_t outputMutex;
extern uint32_t globalTid;
extern uint32_t globalPos;
extern uint32_t globalEnd;
extern uint32_t bin;
extern uint32_t outputBin;
extern uint64_t globalnVariantPositions;

/*! @typedef
 @abstract Global configuration information structure
 @field	keepCpG	1: Output CpG metrics
 @field keepCHG	1: Output CHG metrics
 @field	keepCHH	1: Output CHH metrics
 @field	minMapq	Minimum MAPQ value to include a read (-q)
 @field	minPhred	Minimum Phred score to include a base (-p)
 @field	keepDupes	1: Include marked duplicates when calculating metrics
 @field	minDepth	Minimum depth for outputing methylation metrics
 @field keepDiscordant	0: Do not include discordantly aligned reads when calculating metrics
 @field	keepSingleton	0: Do not include singletons when calculating metrics
 @field ignoreFlags     Mask that's logically &ed with and ignored if > 0. Defaults to 0xF00.
 @field requireFlags    Mask that's logically &ed with and ignored if < mask. Defaults to 0, which means ignore.
 @field merge   1: Merge Cs in either a CpG or CHG context into single entries
 @field methylKit       Output in a format compatible with methylKit
 @field minOppositeDepth	Minimum depth covering the opposite strand needed to look for variants
 @field ignoreNH	If set, don't exclude alignments with NH auxiliary tags with values >1 (i.e., marked multimappers).
 @field maxVariantFrac	If the fraction of non-Gs on the opposite strand is greater than this then a position is excluded.
 @field fraction	1: Output should be the methylation fraction only, 0: otherwise
 @field counts	1: Output just the coverage
 @field logit	1: Logit transform the methylation fraction
 @field cytosine_report	1: Output a bismark-like cytosine report
 @field output_fp	Output file pointers (to CpG, CHG, and CHH, respectively)
 @field	reg	A region specified by -r
 @field BAMName	The BAM file name
 @field minConversionEfficiency	The minimum acceptable conversion efficiency in non-CpG positions for read inclusion
 @field fp	Input file pointer (must be a BAM or CRAM file)
 @field	bai	The index for fp
 @field bedName	The BED file name
 @field bed	Pointer to regions specified in a BED file (-l option)
 @field FastaName	The fasta file name
 @field bounds	Trimming bounds
 @field absoluteBounds	Absolute trimming bounds
 @field nThreads	Number of threads in use.
 @field chunkSize	The number of bases processed by each thread at a time (can be adjusted a bit to ensure CpGs/CHGs aren't split between processors)
*/
typedef struct {
    int keepCpG, keepCHG, keepCHH;
    int minMapq, minPhred, keepDupes, minDepth;
    int keepDiscordant, keepSingleton, ignoreFlags, requireFlags;
    int merge, methylKit, minOppositeDepth;
    int ignoreNH;
    double maxVariantFrac;
    int fraction, counts, logit;
    int cytosine_report;
    FILE **output_fp;
    char *reg;
    char *BAMName;
    float minConversionEfficiency;
    char *BWName;
    char *outBBMName;
    char *BBMName;
    FILE *BBM_ptr;
    char** chromNames;
    uint32_t chromCount;
    uint32_t* chromLengths;
    char filterMappability;
    char outputBB;
    float mappabilityCutoff;
    int minMappableBases;
    char** bw_data;
    htsFile *fp;
    hts_idx_t *bai;
    char *FastaName;
    int bounds[16];
    int absoluteBounds[16];
    int nThreads;
    char noBAM;
    unsigned long chunkSize;
} Config;

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
int getStrand(bam1_t *b);

//Used internally by the pileup-based functions
int updateMetrics(Config *config, const bam_pileup1_t *plp);

//Used internally to parse things like --OT 0,1,2,3
void parseBounds(char *s2, int *vals, int mult);

// Read-pair overlap handling functions
void *initOlapHash();
void destroyOlapHash(void *ohash);

inline int isCpG(char *seq, int pos, int seqlen) {
    if(pos >= seqlen) return 0;
    if(*(seq+pos) == 'C' || *(seq+pos) == 'c') {
        if(pos+1 == seqlen) return 0;
        if(*(seq+pos+1) == 'G' || *(seq+pos+1) == 'g') return 1;
        return 0;
    } else if(*(seq+pos) == 'G' || *(seq+pos) == 'g') {
        if(pos == 0) return 0;
        if(*(seq+pos-1) == 'C' || *(seq+pos-1) == 'c') return -1;
        return 0;
    }
    return 0;
}

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
    int strand = getStrand(b)-1;
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

int getStrand(bam1_t *b) {
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
