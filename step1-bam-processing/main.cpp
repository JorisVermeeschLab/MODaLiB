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
#include <unordered_set>
#include <string>

#include "params.h"
#include "read_info.h"
#include "fragment.h"
#include "methylation/mm_ml.h"
#include "reference.h"
#include "queue_map.h"
#include "strand.h"
#include "vntr.h"
#include "base.h"
#include "pairing/paired_reads.h"
#include "pairing/single_read.h"

#if HTS_VERSION < 101100
#error "HTSlib version should be >= 1.11"
#endif

#define NDEBUG
static constexpr bool DEBUG_TEXT_OUTPUT = false;



pthread_mutex_t positionMutex;
pthread_mutex_t outputMutex;
uint32_t globalTid = 0;
uint32_t globalPos = 0;
uint32_t globalEnd = 0;
uint32_t bin = 0;
uint32_t outputBin = 0;


void write_fragment(kstring_t *os, Fragment& fragment) {

    if (DEBUG_TEXT_OUTPUT) {
        char str[10000];
        snprintf(str, 10000, "%" PRIu8 "\t%" PRIu32 "\t%" PRIu8 "\t%" PRIu32 "\t%" PRIu16" \t%" PRIu16 "\t%" PRIu16 "\t%" PRIu16 "\t%" PRIu16 "\t%" PRIu16 "\n",
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


void trim_alignment(bam1_t *b, const Parameters& config) {

    int strand = get_strand(b);
    const int* bounds = config.get_bounds(strand);

    // Determine the left and right bounds
    int lb, rb;
    if (b->core.flag & BAM_FREAD1) {
        lb = bounds[0];
        rb = bounds[1];
    } else {
        lb = bounds[2];
        rb = bounds[3];
    }

    // Make sure the leftmost bound does not exceed the size of the read
    lb = (lb < b->core.l_qseq) ? lb : b->core.l_qseq;

    // Trim left part of the read
    uint8_t* qual = bam_get_qual(b);
    if (lb) for (int i = 0; i < lb; i++) qual[i] = 0;

    // Trim right part of the read
    if (rb) for (int i = rb; i < b->core.l_qseq; i++) qual[i] = 0;
}


void process_read(const bam1_t *b, const AlignmentRefSequence& ref, ReadInfo& read_info, uint8_t min_phred, bool with_c_conversion) {

    uint32_t query_pos = 0;  // 0-based relative coordinate
    uint32_t ref_pos = b->core.pos;  // 0-based coordinate (relative to some arbitrary start pos)
    uint32_t* CIGAR = bam_get_cigar(b);
    uint8_t* read_seq = bam_get_seq(b);
    uint8_t* read_qual = bam_get_qual(b);
    int strand = get_strand(b);

    for (int cigar_pos = 0; cigar_pos < b->core.n_cigar; cigar_pos++) {

        int cigar_op_type = bam_cigar_op(CIGAR[cigar_pos]);

        for (int kk = 0; kk < bam_cigar_oplen(CIGAR[cigar_pos]); kk++) {

            bool consume_query = true;
            bool consume_ref = true;

            BaseInfo base_info;

            // Get base from the reference sequence
            char ref_base = ref.get_base(ref_pos);

            if ((cigar_op_type == BAM_CMATCH) || (cigar_op_type == BAM_CEQUAL) ||
                (cigar_op_type == BAM_CDIFF)) {  // M, = or X (alignment match, but can be sequence match or mismatch)

                // Get base from query sequence
                char base = seq_nt16_str[bam_seqi(read_seq, query_pos)];
                base_info.c = base;
                base_info.phred_score = read_qual[query_pos];

                // Check that bases indeed match or mismatch
                if (cigar_op_type == BAM_CEQUAL) assert(base == ref_base);
                else if (cigar_op_type == BAM_CDIFF) assert(base != ref_base);

                // Check if CpG site
                int direction = ref.is_cpg(ref_pos, with_c_conversion);
                if (direction != 0) base_info.is_cpg_site = 1;

                if (with_c_conversion) {  // If unmethylated cytosines have been converted (either chemically or enzymatically)

                    // Check methylation status if CpG site, and check the presence of a putative mutation otherwise.
                    // Note: it is generally admitted that CpG site are symmetrically (un)methylated on both strands.
                    if ((direction == 1) && (strand & 1) == 1) {  // OT or CTOT case
                        if (base == 'C') {
                            base_info.methylated = 1;
                        } else if (base == 'T') {
                            base_info.methylated = 0;
                            base_info.c = 'C';
                            assert(ref_base == 'C');
                        }
                    } else if ((direction == -1) && (strand & 1) == 0) {  // OB or CTOB
                        if (base == 'G') {
                            base_info.methylated = 1;
                        } else if (base == 'A') {
                            base_info.methylated = 0;
                            base_info.c = 'G';
                            assert(ref_base == 'G');
                        }
                    } else {  // Not a CpG site
                        if (base != ref_base) {  // base mismatch

                            // Here, the read is assumed to have been aligned to a 3-letter genome,
                            // like in bwa-meth, BSmooth, Bison or Bismark.
                            // If T on OT/CTOT, could be either a C->U->T conversion or a C->T mutation.
                            // If A on OB/CTOB, could be either a A->U->T->A conversion or a G->A mutation.

                            // If bases mismatch, consider this as a putative mutation
                            // unless it's a putative C->T conversion caused by bisulfite
                            // or enzymatic conversion of unmethylated cytosines
                            bool putative_mutation = true;
                            if (((strand & 1) == 1) && (base == 'T') && (ref_base == 'C')) {
                                putative_mutation = false;
                            } else if (((strand & 1) == 0) && (base == 'A') && (ref_base == 'G')) {
                                putative_mutation = false;
                            }
                            if (putative_mutation) {  // Not a converted cytosine -> maybe SNP (or sequencing error)
                                base_info.putative_snp = 1;
                            } else {  // Ambiguous case: can be either a mutation or a converted cytosine
                                // IMPORTANT NOTE: Because cytosines are orders of magnitude more frequent
                                // than SNPs, we consider all mismatches as being caused by converted cytosines.
                                // TODO: maybe adjust Phred score to account for the probability of a SNP?
                                base_info.methylated = 0;
                                base_info.c = ref_base;  // Convert back the converted cytosine
                            }
                        } else {  // Base match
                            // Do nothing
                        }
                    }
                } else {  // No conversion of unmethylated cytosines
                    // In the absence of bisulfite/enzymatic conversion, any base mismatch
                    // is treated as a putative point mutation (or sequencing error)
                    if (base != ref_base) base_info.putative_snp = 1;
                }
            } else if (cigar_op_type == BAM_CINS) {  // I (insertion to the reference)
                char base = seq_nt16_str[bam_seqi(read_seq, query_pos)];
                base_info.c = base;
                base_info.phred_score = read_qual[query_pos];
                base_info.inserted = 1;
                if (with_c_conversion) {
                    // No reference -> no way to disentangle mutations from converted cytosines
                    // -> uninformative base
                    base_info.uninformative = 1;
                }
                consume_ref = false;
            } else if (cigar_op_type == BAM_CDEL) {  // D (deletion from the reference)
                base_info.deleted = 1;
                consume_query = false;
            } else if (cigar_op_type == BAM_CREF_SKIP) {  // N (skipped region from the reference)
                base_info.uninformative = 1;
                consume_query = false;
            } else if ((cigar_op_type == BAM_CSOFT_CLIP) || (cigar_op_type == BAM_CBACK)) {  // S (soft clipping) or B ("backwards" skip operator)
                // Soft-clipped base: base that is not part of the alignment
                base_info.uninformative = 1;
                base_info.soft_clipped = 1;
                char base = seq_nt16_str[bam_seqi(read_seq, query_pos)];
                base_info.c = base;
                base_info.phred_score = read_qual[query_pos];
                consume_ref = false;
            } else {  // Ignore H (hard clipping) and P (padding) operations
                // Hard-clipped base: base that is not part of the alignment and
                // is not present in the sequence at all.
                consume_query = false;
                consume_ref = false;
            }

            // Go to next base
            if (consume_ref) {
                ref_pos++;
                read_info.add_ref_base(ref_base);
            }
            if (consume_query) {
                query_pos++;
            }

            // If base is unknown, set Phred score to 0
            if (base_info.c == 'N') base_info.uninformative = 1;

            // If quality is too low, mark base as uninformative
            if (base_info.phred_score < min_phred) base_info.uninformative = 1;

            read_info.add_base(ref_pos, base_info);
        }
    }

    assert(query_pos == b->core.l_qseq);
}

void* process_partition(void* config_) {
    Parameters* config = static_cast<Parameters*>(config_);
    uint32_t localBin = 0, localPos = 0, localEnd = 0, localTid = 0, localPos2 = 0;
    bam1_t* b = bam_init1();

    kstring_t* os = (kstring_t*) calloc(1, sizeof(kstring_t));
    if (!os) {
        fprintf(stderr, "Couldn't allocate space the for kstring_t structure in perRead!\n");
        return nullptr;
    }

    // Load reference genome
    faidx_t* fai;
    if ((fai = fai_load(config->fasta_filename)) == nullptr) {
        fprintf(stderr, "Couldn't open the index for %s!\n", config->fasta_filename);
        return nullptr;
    }

    // Open BAM file
    htsFile* fp;
    if ((fp = hts_open(config->bam_filename, "rb")) == nullptr) {
        fprintf(stderr, "Couldn't open BAM file %s\n", config->bam_filename);
        return nullptr;
    }

    // Open BAM index file
    hts_idx_t* bai;
    if ((bai = sam_index_load(fp, config->bam_filename)) == nullptr) {
        fprintf(stderr, "Couldn't open BAM index %s\n", config->bam_filename);
        return nullptr;
    }
    bam_hdr_t* hdr = sam_hdr_read(fp);

    QueueMap<std::string, ReadInfo> paired_reads(100000);
    std::unordered_set<std::string> blacklisted_pairs;

    while (true) {
        //Lock and unlock the mutex, so we can get/update the tid/position
        pthread_mutex_lock(&positionMutex);
        localBin = bin++;
        localTid = globalTid;
        localPos = globalPos;
        localEnd = localPos + config->partition_size;
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
        if (localTid < hdr->n_targets && globalTid != std::numeric_limits<uint32_t>::max()) {
            if (globalPos >= hdr->target_len[localTid]) {
                localEnd = hdr->target_len[localTid];
                globalTid++;
                globalPos = 0;
            }
        }
        pthread_mutex_unlock(&positionMutex);

        localPos2 = 0;
        if (localPos > 1) localPos2 = localPos - 2;

        // Break out of the loop if finished
        if (localTid >= hdr->n_targets) break;
        if (globalEnd && localPos >= globalEnd) break;

        hts_itr_t* iter = sam_itr_queryi(bai, localTid, localPos, localEnd);

        //Fetch the sequence and then iterate over the reads. There's a 10kb buffer on the end
        int seq_len;
        char* seq = faidx_fetch_seq(fai, hdr->target_name[localTid], localPos2, localEnd+10000, &seq_len);
        while (sam_itr_next(fp, iter, b) >= 0) {
            if (b->core.pos < localPos) continue;
            if (b->core.pos >= localEnd) break;

            // Create read
            std::string id = std::string(bam_get_qname(b));
            ReadInfo read(
                id, hdr->target_name[b->core.tid],
                b->core.pos + 1,  // Change 0-based to 1-based coordinate
                abs(b->core.l_qseq), abs(b->core.isize),
                get_strand(b),
                b->core.flag
            );

            // If secondary alignment, or mapping quality is insufficient, discard read.
            if ((!config->single_read) or (blacklisted_pairs.find(read.id) != blacklisted_pairs.end())) {
                // The other read in the pair has been discarded due to low mapping quality
                // -> Discard the current one as well.
                blacklisted_pairs.erase(read.id);
                continue;
            } else if ((b->core.qual < config->min_mapq) or (!read.is_worth_analyzing())) {
                if (!config->single_read) {
                    // The other read (upcoming) should be discarded as well
                    blacklisted_pairs.insert(read.id);
                }
                continue;
            }

            // Use inclusion bounds to trim alignments
            trim_alignment(b, *config);

            // Extract read-level statistics from read
            AlignmentRefSequence ref(seq, localPos2, seq_len);
            process_read(
                    b, ref, read, config->min_phred,
                    !config->use_mm_ml  // IMPORTANT: When using MM/ML flags to call methylation, we assume that Cs have not been converted
            );

            // Call methylation using MM and ML tags (if required)
            if (config->use_mm_ml) methylation_from_aux(b, read, config->ml_cutoff);
            // TODO: call methylation from aligner-specific optional fields such as Bismark (XM tag), if using them

            // Combine read pairs into fragments
            if (config->single_read) {  // Single-read
                Fragment fragment = compute_fragment_statistics(read);
                if (fragment.is_valid()) write_fragment(os, fragment);
            } else {  // Paired-end
                bool found_pair = false;
                if (paired_reads.contains(read.id)) {
                    ReadInfo& other_read = paired_reads.get(read.id);
                    if (read.paired_with(other_read)) {
                        Fragment fragment = compute_fragment_statistics(
                                read, other_read, ref, config->min_phred);
                        if (fragment.is_valid()) write_fragment(os, fragment);
                        found_pair = true;
                        paired_reads.remove(read.id);
                    }
                }
                if (!found_pair) paired_reads.add(read.id, read);
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
"\n"
"Usage: BAMExtractor [OPTIONS] <ref.fa> <input>\n"
"\n"
"Positional arguments:\n"
"  ref.fa    Reference genome in fasta format. This must be indexed with\n"
"            samtools faidx.\n"
"  input     An input BAM or CRAM file. This MUST be sorted and should be indexed.\n"
"\n"
"Optional arguments:\n"
" -q INT     Minimum MAPQ threshold to include an alignment (default 10)\n"
" -p INT     Minimum Phred threshold to include a base (default 5). This must "
"            be >0.\n"
" -o STR     Output file name [stdout]\n"
" -nthreads INT    The number of threads to use, the default is 1\n"
" --partition_size INT  The size of the genome processed by a single thread at a time.\n"
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
" --OB INT,INT,INT,INT As with --OT, but for the original bottom.\n"
" --CTOT INT,INT,INT,INT As with --OT, but for complementary to the original top.\n"
" --CTOB INT,INT,INT,INT As with --OT, but for complementary to the original bottom.\n"
"\n"
"Note that this program will produce incorrect values for alignments spanning\n"
"more than 10kb.\n");
}

int run_convertor(int argc, char *argv[]) {
    Parameters config;
    FILE *ofile = stdout;
    faidx_t *fai;

    static struct option lopts[] = {
            {"help",    0, nullptr, 'h'},
            {"version", 0, nullptr, 'v'},
            {"single_read",   0, nullptr, 20},
            {"useMMML",   0, nullptr, 21},
            {"mlCutoff",   1, nullptr, 22},
            {"partition_size",    1, nullptr,  19},
            {"OT",           1, nullptr,   7},
            {"OB",           1, nullptr,   8},
            {"CTOT",         1, nullptr,   9},
            {"CTOB",         1, nullptr,  10},
            {0,         0, nullptr,   0}
    };
    char c;
    while((c = getopt_long(argc, argv, "hvq:p:o:@:r:l:F:R:", lopts, nullptr)) >= 0) {
        switch(c) {
            case 'h':
                cli_help();
                return 0;
            case 'o':
                if((ofile = fopen(optarg, "wb")) == nullptr) {
                    fprintf(stderr, "Couldn't open %s for writing\n", optarg);
                    return 2;
                }
                break;
            case 'q':
                config.min_mapq = atoi(optarg);
                break;
            case 'p':
                config.min_phred = atoi(optarg);
                break;
            case '@':
                config.n_threads = atoi(optarg);
                break;
            case 20:
                config.single_read = 1;
                break;
            case 21:
                config.use_mm_ml = 1;
                break;
            case 22:
                config.ml_cutoff = static_cast<uint8_t>(atoi(optarg));
                break;
            case 7:
                config.parse_ot_bounds(optarg);
                break;
            case 8:
                config.parse_ob_bounds(optarg);
                break;
            case 9:
                config.parse_ctot_bounds(optarg);
                break;
            case 10:
                config.parse_ctob_bounds(optarg);
                break;
            case 19:
                config.partition_size = strtoul(optarg, nullptr, 10);
                if(config.partition_size < 1) {
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
    if (config.min_phred < 1) {
        fprintf(stderr, "-p %i is invalid. resetting to 1, which is the lowest possible value.\n", config.min_phred);
        config.min_phred = 1;
    }
    if (config.min_mapq < 0) {
        fprintf(stderr, "-q %i is invalid. Resetting to 0, which is the lowest possible value.\n", config.min_mapq);
        config.min_mapq = 0;
    }

    if ((fai = fai_load(argv[optind])) == nullptr) {
        fprintf(stderr, "Couldn't open the index for %s!\n", argv[optind]);
        cli_help();
        return -2;
    }
    config.fasta_filename = argv[optind];
    config.bam_filename = argv[optind + 1];
    if ((config.fp = hts_open(argv[optind+1], "rb")) == nullptr) {
        fprintf(stderr, "Couldn't open %s for reading!\n", argv[optind+1]);
        return -4;
    }
    if ((config.bai = sam_index_load(config.fp, argv[optind+1])) == nullptr) {
        fprintf(stderr, "Couldn't load the index for %s, will attempt to build it.\n", argv[optind+1]);
        if (bam_index_build(argv[optind+1], 0) < 0) {
            fprintf(stderr, "Couldn't build the index for %s! File corrupted?\n", argv[optind+1]);
            return -5;
        }
        if ((config.bai = sam_index_load(config.fp, argv[optind+1])) == nullptr) {
            fprintf(stderr, "Still couldn't load the index, quiting.\n");
            return -5;
        }
    }

    // Output file
    config.output_fp = (FILE**) malloc(sizeof(FILE *));
    config.output_fp[0] = ofile;

    // Spin up the threads
    pthread_mutex_init(&positionMutex, nullptr);
    pthread_t *threads = (pthread_t*) calloc(config.n_threads, sizeof(pthread_t));
    for (int i = 0; i < config.n_threads; i++) pthread_create(threads + i, nullptr, &process_partition, &config);
    for (int i = 0; i < config.n_threads; i++) pthread_join(threads[i], nullptr);
    free(threads);
    pthread_mutex_destroy(&positionMutex);

    // Close files and free memory
    hts_close(config.fp);
    if (ofile != stdout) fclose(ofile);
    free(config.output_fp);
    fai_destroy(fai);
    hts_idx_destroy(config.bai);

    return 0;
}

int main(int argc, char *argv[]) {
    if ((argc == 1) || (strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "--help") == 0)) {
        cli_help();
        return 0;
    } else {
        return run_convertor(argc, argv);
    }
}
