// IMPORTANT NOTE: this file contains source code from multiple files taken from the MethylDackel software

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

#include <getopt.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <pthread.h>
#include <limits>
#include <string>
#include <vector>

#include "MethylDackel.h"
#include "read_info.h"
#include "fragment.h"
#include "methylation.h"
#include "utils.h"


#define NDEBUG
static constexpr bool DEBUG_TEXT_OUTPUT = true;



pthread_mutex_t positionMutex;
pthread_mutex_t outputMutex;
uint32_t globalTid = 0;
uint32_t globalPos = 0;
uint32_t globalEnd = 0;
uint32_t bin = 0;
uint32_t outputBin = 0;


void addFragment(kstring_t *os, Fragment& fragment) {

    if (fragment.start_chr_id == fragment.end_chr_id) {
        // Quality checks
        if (!fragment.mapped_in_proper_pair()) return;
        if (fragment.failed_quality_checks()) return;
        if (fragment.is_duplicate()) return;
    } else {
        // Do nothing: keep fragments that may potentially highlight translocations
    }

    // Ensure known chromosome
    if ((fragment.start_chr_id == 0) || (fragment.end_chr_id == 0)) return;

    // Ensure fragment is not too long
    // if (read.insert_size > 1000) return;

    if (DEBUG_TEXT_OUTPUT) {
        char str[10000];
        snprintf(str, 10000, "COMBINED: %" PRIu8 "\t%" PRIu32 "\t%" PRIu8 "\t%" PRIu32 "\t%" PRIu16" \t%" PRIu16 "\t%" PRIu16 "\t%" PRIu16 "\t%" PRIu16 "\t%" PRIu16 "\n",
            fragment.start_chr_id,
            fragment.start_pos,
            fragment.end_chr_id,
            fragment.end_pos,
            fragment.length,
            fragment.n_methylated,
            fragment.n_unmethylated,
            fragment.snp_count,
            fragment.gc_count,
            fragment.s_count
        );
        kputs(str, os);
    } else {
        char buffer[FRAGMENT_SIZE_IN_BYTES];
        fragment.serialize(buffer);
        kputsn(buffer, sizeof(buffer), os);
    }
}


void addRead(kstring_t *os, bam1_t *b, bam_hdr_t *hdr, ReadInfo& read) {

    // Do some checks
    if (read.failed_quality_checks()) return;
    if (read.is_duplicate()) return;
    if (!read.mapped_in_proper_pair()) return;

    // Ensure known chromosome
    if (read.chr_id == 0) return;

    char str[10000];

    snprintf(str, 10000, "%" PRId32 "\t%" PRIu16 "\t%" PRId32 "\t%" PRIu32 "\t%" PRIu32 "\t%" PRIu32 "\t%" PRIu32 "\t%" PRId16 "\t%" PRIu16 "\t%" PRIu16 "\t%" PRIu16 "\n",
        read.length,
        read.chr_id,
        read.start,
        read.n_methylated,
        read.n_unmethylated,
        read.gc_count,
        read.s_count,
        read.insert_size,
        read.sam_flag,
        read.inserted,
        read.deleted
    );

    kputs(str, os);
}


// TODO: replace seq, sequenceStart and seqLen by a custom object ("AlignedRefSequence"?)
void compute_read_statistics(bam1_t *b, char *seq, uint32_t sequenceStart, int seqLen, ReadInfo& read_info, int min_phred) {

    uint32_t query_pos = 0;
    uint32_t ref_pos = b->core.pos;
    uint32_t* CIGAR = bam_get_cigar(b);
    uint8_t* read_seq = bam_get_seq(b);
    uint8_t* read_qual = bam_get_qual(b);
    int strand = get_strand(b);

    while (int cigar_pos = 0; cigar_pos < b->core.n_cigar; cigar_pos++) {

        int cigar_op_type = bam_cigar_op(CIGAR[cigar_pos]);

        if ((cigar_op_type == BAM_CMATCH) || (cigar_op_type == BAM_CEQUAL) || (cigar_op_type == BAM_CDIFF)) {  // M, = or X (alignment match, but can be sequence match or mismatch)

            for (int kk = 0; kk < bam_cigar_oplen(CIGAR[cigar_pos]); kk++) {
                if (read_qual[query_pos] >= min_phred) {  // Consider only high-quality base calls

                    // Get bases from the reference and query sequences
                    int base = bam_seqi(read_seq, query_pos);
                    int ref_base = *(seq + ref_pos - sequenceStart);

                    // Check that bases indeed match or mismatch
                    if (cigar_op_type == BAM_CEQUAL) assert(base == ref_base);
                    else if (cigar_op_type == BAM_CDIFF) assert(base != ref_base);

                    // Update read-level statistics with the new base
                    bool is_snp = ((base != ref_base) && (static_cast<char>(ref_base) != 'C'));
                    // TODO: count non_snp bases
                    read_info.add_base(base, query_pos, strand, is_snp);

                    // Check if CpG site
                    int direction = isCpG(seq, ref_pos - sequenceStart, seqLen);

                    // Check methylation status (if CpG)
                    if (direction == 1 && (strand & 1) == 1) { // C & OT/CTOT
                        if (base == 2) read_info.inc_methyl_count(query_pos);  //C
                        else if (base == 8) read_info.inc_unmethyl_count(query_pos);  //T
                        read_info.cpgs.push_back(query_pos);
                    } else if (direction == -1 && (strand & 1) == 0) { // G & OB/CTOB
                        if (base == 4) read_info.inc_methyl_count(query_pos);  //G
                        else if (base == 1) read_info.inc_unmethyl_count(query_pos);  //A
                        read_info.cpgs.push_back(query_pos);
                    }
                }
                ref_pos++;
                query_pos++;
            }
        } else if (cigar_op_type == BAM_CINS) {  // I (insertion to the reference)

            // Update read statistics with the inserted bases
            int inserted = bam_cigar_oplen(CIGAR[cigar_pos]);
            for (int i = 0; i < inserted; i++) {
                int base = bam_seqi(read_seq, query_pos);
                read_info.add_base(base, query_pos, strand);
                query_pos++;
            }

            read_info.inserted += inserted;
        } else if (cigar_op_type == BAM_CDEL) {  // D (deletion from the reference)
            int n_deleted = bam_cigar_oplen(CIGAR[cigar_pos]);
            read_info.deleted += n_deleted;
            ref_pos += n_deleted;
        } else if (cigar_op_type == BAM_CREF_SKIP) {  // N (skipped region from the reference)
            int n_skipped = bam_cigar_oplen(CIGAR[cigar_pos]);
            ref_pos += n_skipped;
        } else if ((cigar_op_type == BAM_CSOFT_CLIP) || (cigar_op_type == BAM_CBACK)) {  // S (soft clipping) or B ("backwards" skip operator)
            int n_clipped = bam_cigar_oplen(CIGAR[cigar_pos]);
            query_pos += n_clipped;
            cigar_n_bases_processed += n_clipped;
        } else {  // Ignore H and P operations
            // Do nothing
        }
    }

    assert(query_pos == b->core.l_qseq);
}

void* perReadMetrics(void *foo) {
    Config *config = (Config*) foo;
    bam_hdr_t *hdr;
    int seqlen;
    uint32_t localBin = 0, localPos = 0, localEnd = 0, localTid = 0, localPos2 = 0;
    char *seq = NULL;
    kstring_t *os = NULL;
    faidx_t *fai;
    hts_idx_t *bai;
    htsFile *fp;
    bam1_t *b = bam_init1();
    hts_itr_t *iter;

    os = (kstring_t*) calloc(1, sizeof(kstring_t));
    if(!os) {
        fprintf(stderr, "Couldn't allocate space the for kstring_t structure in perRead!\n");
        return nullptr;
    }

    //Open the files
    if((fai = fai_load(config->FastaName)) == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", config->FastaName);
        return nullptr;
    }
    if((fp = hts_open(config->BAMName, "rb")) == NULL) {
        fprintf(stderr, "Couldn't open %s for reading!\n", config->BAMName);
        return nullptr;
    }
    if((bai = sam_index_load(fp, config->BAMName)) == NULL) {
        fprintf(stderr, "Couldn't load the index for %s\n", config->BAMName);
        return nullptr;
    }
    hdr = sam_hdr_read(fp);

    std::vector<ReadInfo> paired_reads;

    while (true) {
        //Lock and unlock the mutex, so we can get/update the tid/position
        pthread_mutex_lock(&positionMutex);
        localBin = bin++;
        localTid = globalTid;
        localPos = globalPos;
        localEnd = localPos + config->chunkSize;
        if(localTid >= hdr->n_targets) {
            pthread_mutex_unlock(&positionMutex);
            break;
        }
        if (globalEnd && localEnd > globalEnd) localEnd = globalEnd;
        globalPos = localEnd;
        if (globalEnd > 0 && globalPos >= globalEnd) {
            //If we've specified a region, then break once we're outside of it
            globalTid = std::numeric_limits<uint32_t>::max();
        }
        if(localTid < hdr->n_targets && globalTid != std::numeric_limits<uint32_t>::max()) {
            if(globalPos >= hdr->target_len[localTid]) {
                localEnd = hdr->target_len[localTid];
                globalTid++;
                globalPos = 0;
            }
        }
        pthread_mutex_unlock(&positionMutex);

        localPos2 = 0;
        if(localPos > 1) {
            localPos2 = localPos - 2;
        }

        //Break out of the loop if finished
        if(localTid >= hdr->n_targets) break; //Finish looping
        if(globalEnd && localPos >= globalEnd) break;
        iter = sam_itr_queryi(bai, localTid, localPos, localEnd);

        //Fetch the sequence and then iterate over the reads. There's a 10kb buffer on the end
        seq = faidx_fetch_seq(fai, hdr->target_name[localTid], localPos2, localEnd+10000, &seqlen);
        while (sam_itr_next(fp, iter, b) >= 0) {
            if (b->core.pos < localPos) continue;
            if (b->core.pos >= localEnd) break;

            // Create read
            std::string id = std::string(bam_get_qname(b));
            ReadInfo read(
                id, hdr->target_name[b->core.tid], b->core.pos,
                abs(b->core.l_qseq), abs(b->core.isize), b->core.flag
            );

            if (config->requireFlags && (config->requireFlags & b->core.flag) != config->requireFlags) continue;
            if (config->ignoreFlags && (config->ignoreFlags & b->core.flag) != 0) continue;
            if (b->core.qual < config->minMapq) continue;

            // Use inclusion bounds to trim alignments
            b = trimAlignmentSimple(b, config->bounds);

            // Extract read-level statistics from read
            compute_read_statistics(b, seq, localPos2, seqlen, read, config->minPhred);

            // Compute methylation using MM and ML tags (if required)
            if (config->use_mm_ml) {
                std::pair<int, int> meth = methylation_from_aux(b, config->ml_cutoff);
                read.n_methylated = meth.first;
                read.n_methylated_nov = meth.first;
                read.n_unmethylated = meth.second;
                read.n_unmethylated_nov = meth.second;
            }

            // Combine read pairs into fragments
            if (config->singleRead) {  // Single-read
                Fragment fragment(read);
                addFragment(os, fragment);
            } else {  // Paired-end
                bool found_pair = false;
                auto pr_iter = paired_reads.begin();
                while (pr_iter != paired_reads.end()) {
                    if (read.paired_with(*pr_iter)) {
                        paired_reads.erase(pr_iter);
                        ReadInfo r1, r2;
                        if (pr_iter->start < read.start) {
                            r1 = *pr_iter;
                            r2 = read;
                        } else {
                            r1 = read;
                            r2 = *pr_iter;
                        }
                        Fragment fragment(r1, r2);
                        if (DEBUG_TEXT_OUTPUT) {
                            addRead(os, b, hdr, r1);
                            addRead(os, b, hdr, r2);
                        }
                        addFragment(os, fragment);
                        found_pair = true;
                        break;
                    }
                    ++pr_iter;
                }
                if (!found_pair) {
                    paired_reads.push_back(read);
                    if (paired_reads.size() > 40) {
                        paired_reads.erase(paired_reads.begin());
                    }
                }
            }

        }
        sam_itr_destroy(iter);
        free(seq);

        // Write output
        while (true) {
            pthread_mutex_lock(&outputMutex);
            if (outputBin != localBin) {
                pthread_mutex_unlock(&outputMutex);
                continue;
            }
            if (os->l) fwrite(os->s, 1, os->l, config->output_fp[0]);
            os->l = 0;
            outputBin++;
            pthread_mutex_unlock(&outputMutex);
            break;
        }
    }

    if (os->s) free(os->s);
    free(os);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    fai_destroy(fai);
    hts_close(fp);
    hts_idx_destroy(bai);
    return nullptr;
}

void cli_help() {
    fprintf(stderr,
"\nUsage: BAMExtractor [OPTIONS] <ref.fa> <input>\n"
"\n"
"Arguments:\n"
"  ref.fa    Reference genome in fasta format. This must be indexed with\n"
"            samtools faidx\n"
"  input     An input BAM or CRAM file. This MUST be sorted and should be indexed.\n"
"\nOptions:\n"
" -q INT     Minimum MAPQ threshold to include an alignment (default 10)\n"
" -p INT     Minimum Phred threshold to include a base (default 5). This must "
"            be >0.\n"
" -r STR     Region string in which to extract methylation\n"
" -o STR     Output file name [stdout]\n"
" -F         By default, all reads are output. If you would like to\n"
"            ignore certain classes of reads then simply give a value for their\n"
"            flags here. Note that an alignment will be logically anded with this\n"
"            flag, so a single bit overlap will lead to exclusion. The default\n"
"            for this is 0.\n"
" -R         Require each alignment to have all bits in this value\n"
"            present, or else the alignment is ignored. This is equivalent to the\n"
"            -f option in samtools. The default is 0, which includes all\n"
"            alignments.\n"
" -nthreads INT    The number of threads to use, the default is 1\n"
" --chunkSize INT  The size of the genome processed by a single thread at a time.\n"
"            The default is 1000000 bases. This value MUST be at least 1.\n"
" --OT INT,INT,INT,INT Inclusion bounds for methylation calls from reads/pairs\n"
"                  originating from the original top strand. Suggested values can\n"
"                  be obtained from the MBias program. Each integer represents a\n"
"                  1-based position on a read. For example --OT A,B,C,D\n"
"                  translates to, \"Include calls at positions from A through B\n"
"                  on read #1 and C through D on read #2\". If a 0 is used a any\n"
"                  position then that is translated to mean start/end of the\n"
"                  alignment, as appropriate. For example, --OT 5,0,0,0 would\n"
"                  include all but the first 4 bases on read #1. Users are\n"
"                  strongly advised to consult a methylation bias plot, for\n"
"                  example by using the MBias program.\n"
" --OB INT,INT,INT,INT\n"
" --CTOT INT,INT,INT,INT\n"
" --CTOB INT,INT,INT,INT As with --OT, but for the original bottom, complementary\n"
"                  to the original top, and complementary to the original bottom\n"
"                  strands, respectively.\n"
"\n"
"Note that this program will produce incorrect values for alignments spanning\n"
"more than 10kb.\n");
}

int run_convertor(int argc, char *argv[]) {
    Config config;
    FILE *ofile = stdout;
    int i;
    faidx_t *fai;
    bam_hdr_t *hdr = NULL;
    char c;

    static struct option lopts[] = {
            {"help",    0, NULL, 'h'},
            {"version", 0, NULL, 'v'},
            {"singleRead",   0, NULL, 20},
            {"useMMML",   0, NULL, 21},
            {"mlCutoff",   1, NULL, 22},
            {"chunkSize",    1, NULL,  19},
            {"ignoreFlags",  1, NULL, 'F'},
            {"requireFlags", 1, NULL, 'R'},
            {"OT",           1, NULL,   7},
            {"OB",           1, NULL,   8},
            {"CTOT",         1, NULL,   9},
            {"CTOB",         1, NULL,  10},
            {0,         0, NULL,   0}
    };
    while((c = getopt_long(argc, argv, "hvq:p:o:@:r:l:F:R:", lopts, NULL)) >= 0) {
        switch(c) {
            case 'h':
                cli_help();
                return 0;
            case 'o':
                if((ofile = fopen(optarg, "wb")) == NULL) {
                    fprintf(stderr, "Couldn't open %s for writing\n", optarg);
                    return 2;
                }
                break;
            case 'q':
                config.minMapq = atoi(optarg);
                break;
            case 'p':
                config.minPhred = atoi(optarg);
                break;
            case '@':
                config.nThreads = atoi(optarg);
                break;
            case 'r':
                config.reg = optarg;
                break;
            case 'F':
                config.ignoreFlags = atoi(optarg);
                break;
            case 'R':
                config.requireFlags = atoi(optarg);
                break;
            case 20:
                config.singleRead = 1;
                break;
            case 21:
                config.use_mm_ml = 1;
                break;
            case 22:
                config.ml_cutoff = static_cast<uint8_t>(atoi(optarg));
                break;
            case 7:
                parseBounds(optarg, config.bounds, 0);
                break;
            case 8:
                parseBounds(optarg, config.bounds, 1);
                break;
            case 9:
                parseBounds(optarg, config.bounds, 2);
                break;
            case 10:
                parseBounds(optarg, config.bounds, 3);
                break;
            case 19:
                config.chunkSize = strtoul(optarg, NULL, 10);
                if(config.chunkSize < 1) {
                    fprintf(stderr, "Error: The chunk size must be at least 1!\n");
                    return 1;
                }
                break;
            default:
                fprintf(stderr, "Invalid option '%c'\n", c);
                cli_help();
                return 1;
        }
    }

    if (argc == 1) {
        cli_help();
        return 0;
    }
    if (argc - optind != 2) {
        fprintf(stderr, "You must supply a reference genome in fasta format and a BAM or CRAM file\n");
        cli_help();
        return -1;
    }

    // Are the options reasonable?
    if (config.minPhred < 1) {
        fprintf(stderr, "-p %i is invalid. resetting to 1, which is the lowest possible value.\n", config.minPhred);
        config.minPhred = 1;
    }
    if (config.minMapq < 0) {
        fprintf(stderr, "-q %i is invalid. Resetting to 0, which is the lowest possible value.\n", config.minMapq);
        config.minMapq = 0;
    }

    if ((fai = fai_load(argv[optind])) == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", argv[optind]);
        cli_help();
        return -2;
    }
    config.FastaName = argv[optind];
    config.BAMName = argv[optind+1];
    if ((config.fp = hts_open(argv[optind+1], "rb")) == NULL) {
        fprintf(stderr, "Couldn't open %s for reading!\n", argv[optind+1]);
        return -4;
    }
    if ((config.bai = sam_index_load(config.fp, argv[optind+1])) == NULL) {
        fprintf(stderr, "Couldn't load the index for %s, will attempt to build it.\n", argv[optind+1]);
        if (bam_index_build(argv[optind+1], 0) < 0) {
            fprintf(stderr, "Couldn't build the index for %s! File corrupted?\n", argv[optind+1]);
            return -5;
        }
        if ((config.bai = sam_index_load(config.fp, argv[optind+1])) == NULL) {
            fprintf(stderr, "Still couldn't load the index, quiting.\n");
            return -5;
        }
    }

    // Output file
    config.output_fp = (FILE**) malloc(sizeof(FILE *));
    config.output_fp[0] = ofile;

    // Parse the region, if needed
    if (config.reg) {
        const char *foo;
        char *bar = NULL;
        int s, e;
        hdr = sam_hdr_read(config.fp);
        foo = hts_parse_reg(config.reg, &s, &e);
        if(foo == NULL) {
            fprintf(stderr, "Could not parse the specified region!\n");
            return -4;
        }
        bar = (char*) malloc(foo - config.reg + 1);
        if(bar == NULL) {
            fprintf(stderr, "Could not allocate temporary space for parsing the requested region!\n");
            return -5;
        }
        strncpy(bar, config.reg, foo - config.reg);
        bar[foo - config.reg] = 0;
        globalTid = bam_name2id(hdr, bar);
        if(globalTid == (uint32_t) -1) {
            fprintf(stderr, "%s did not match a known chromosome/contig name!\n", config.reg);
            return -6;
        }
        if(s>0) globalPos = s;
        if(e>0) globalEnd = e;
        if(globalEnd > hdr->target_len[globalTid]) globalEnd = hdr->target_len[globalTid];
        free(bar);
    }

    // Spin up the threads
    pthread_mutex_init(&positionMutex, NULL);
    pthread_t *threads = (pthread_t*) calloc(config.nThreads, sizeof(pthread_t));
    for (i=0; i < config.nThreads; i++) pthread_create(threads+i, NULL, &perReadMetrics, &config);
    for (i=0; i < config.nThreads; i++) pthread_join(threads[i], NULL);
    free(threads);
    pthread_mutex_destroy(&positionMutex);

    // Close files and free memory
    hts_close(config.fp);
    if(ofile != stdout) fclose(ofile);
    free(config.output_fp);
    fai_destroy(fai);
    hts_idx_destroy(config.bai);

    return 0;
}

int main(int argc, char *argv[]) {
    if (argc == 1) {
        cli_help();
        return 0;
    } else if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        cli_help();
        return 0;
    } else {
        return run_convertor(argc, argv);
    }
}
