// IMPORTANT NOTE: this file contains source code from multiple files taken from the MethylDackel software

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

#include <getopt.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <cinttypes>
#include <regex>
#include <pthread.h>
#include <limits>
#include <algorithm>
#include <unordered_set>
#include <map>
#include "msi/allele_length.h"
#include <atomic>
#include <cctype>
#include <string>

#include "params.h"
#include "read_info.h"
#include "fragment.h"
#include "methylation/mm_ml.h"
#include "reference.h"
#include "queue_map.h"
#include "strand.h"
#include "base.h"
#include "pairing/paired_reads.h"
#include "pairing/single_read.h"
#include "utils.h"

#if HTS_VERSION < 101100
#error "HTSlib version should be >= 1.11"
#endif

// #define NDEBUG


pthread_mutex_t positionMutex;
pthread_mutex_t outputMutex;
pthread_mutex_t msiMutex = PTHREAD_MUTEX_INITIALIZER;
static std::atomic<uint64_t> msi_reads_checked{0};
static std::atomic<uint64_t> msi_candidates_found{0};
static std::atomic<uint64_t> msi_records_written{0};
uint32_t globalTid = 0;
uint32_t globalPos = 0;
uint32_t globalEnd = 0;
uint32_t bin = 0;
uint32_t outputBin = 0;


void write_fragment(kstring_t *os, const Fragment& fragment) {
    // Save fragment information in binary file
    char buffer[FRAGMENT_SIZE_IN_BYTES];
    fragment.serialize(buffer);
    kputsn(buffer, sizeof(buffer), os);
}

// MSI output: a reference-anchored locus plus the observed sample allele.
// First four columns are preserved for backward compatibility.
// Fractional repeat counts are ratios (bp / motif length), not motif-purity calls.
static void write_msi_record(FILE* fp, uint8_t chr_id,
                             const MicrosatelliteAllele& allele,
                             float weight) {
    if (!fp || allele.start == 0 || allele.motif.empty()) return;

    pthread_mutex_lock(&msiMutex);

    fprintf(fp, "%u\t%u\t%u\t%u\t%s\t%u\t%.3f\t%.3f\t%.3f\n",
            static_cast<unsigned int>(chr_id),
            allele.start - 1,
            allele.end,
            allele.observed_length,
            allele.motif.c_str(),
            allele.reference_length,
            allele.observed_repeat_count(),
            allele.reference_repeat_count(),
            static_cast<double>(weight));

    pthread_mutex_unlock(&msiMutex);
    ++msi_records_written;
}

// Detect reference-defined repeat loci, then measure the observed length
// in this sample's aligned read bases (including indels).
static void detect_microsatellites_in_read(
    FILE* fp,
    const ReadInfo& read,
    uint32_t flank,
    float weight
) {
    if (!fp || !read.chr_id) return;

    ++msi_reads_checked;

    const auto calls = msi_call_alleles(read, flank);
    msi_candidates_found += calls.size();

    for (const auto& call : calls) {
        write_msi_record(fp, read.chr_id, call, weight);
    }
}

static void detect_microsatellites_in_pair(FILE* fp, const ReadInfo& r1,
                                           const ReadInfo& r2, uint32_t flank,
                                           float weight) {
    if (!fp || !r1.chr_id || r1.chr_id != r2.chr_id) return;
    msi_reads_checked += 2;
    const auto calls1 = msi_call_alleles(r1, flank);
    const auto calls2 = msi_call_alleles(r2, flank);
    msi_candidates_found += calls1.size() + calls2.size();

    // A locus seen by both mates is reported once. If their lengths disagree,
    // discard the locus rather than silently choosing an arbitrary allele.
    std::map<std::string, MicrosatelliteAllele> by_locus;
    std::unordered_set<std::string> conflicting;
    auto add = [&](const MicrosatelliteAllele& call) {
        const std::string key = std::to_string(call.start) + ":" +
            std::to_string(call.end) + ":" + call.motif;
        auto result = by_locus.emplace(key, call);
        if (!result.second && result.first->second.observed_length != call.observed_length)
            conflicting.insert(key);
    };
    for (const auto& call : calls1) add(call);
    for (const auto& call : calls2) add(call);
    for (const auto& pair : by_locus) {
        if (conflicting.count(pair.first)) continue;
        const auto& call = pair.second;
        write_msi_record(fp, r1.chr_id, call, weight);
    }
}

void write_methylation(kstring_t *c_os, OrderedQueueMap<genome_size_t, std::tuple<uint8_t, uint16_t, uint16_t, uint8_t>>& queue, const Fragment& fragment) {

    for (const auto& p : fragment.cytosines) {
        const BaseInfo& base = p.second;
        genome_size_t coords = combine_chromosomal_coordinates(fragment.chr_id, p.first);
        if ((base.context == CpG_CONTEXT) && ((fragment.strand & 1) == 0)) coords--;
        uint8_t c_type = as_uint8_t(base.context);

        if (!queue.contains(coords)) {
            queue.add(coords, std::make_tuple(c_type, 0, 0, fragment.strand));
        }
        if (base.methylated) std::get<1>(queue.get(coords))++;
        else std::get<2>(queue.get(coords))++;
    }

    // Save methylation information to file
    while (queue.size_exceeded() && (queue.size() > 0)) {
        auto p = queue.front();
        std::pair<uint8_t, uint32_t> coords = uncombine_chromosomal_coordinates(p.first);
        queue.pop();

        // Serialize methylation information
        int offset = 0;
        char buffer[11];
        offset += Serializer::serialize(buffer + offset, (uint8_t) coords.first);
        offset += Serializer::serialize(buffer + offset, (uint32_t) coords.second);
        offset += Serializer::serialize(buffer + offset, (uint8_t) std::get<0>(p.second));
        offset += Serializer::serialize(buffer + offset, (uint16_t) std::get<1>(p.second));
        offset += Serializer::serialize(buffer + offset, (uint16_t) std::get<2>(p.second));
        offset += Serializer::serialize(buffer + offset, (uint8_t) std::get<3>(p.second));
        kputsn(buffer, sizeof(buffer), c_os);
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


void process_read(const bam1_t *b, const AlignmentRefSequence& ref, ReadInfo& read_info, const Parameters* params) {

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
            bool valid_base = true;

            BaseInfo base_info;

            // Get base from the reference sequence
            char ref_base = ref.get_base(ref_pos);
            BaseInfo ref_base_info = BaseInfo();
            ref_base_info.c = ref_base;

            if ((cigar_op_type == BAM_CMATCH) || (cigar_op_type == BAM_CEQUAL) ||
                (cigar_op_type == BAM_CDIFF)) {  // M, = or X (alignment match, but can be sequence match or mismatch)

                // Get bases from query sequence
                char base = seq_nt16_str[bam_seqi(read_seq, query_pos)];  // Base at current position
                base_info.c = base;
                base_info.phred_score = read_qual[query_pos];
                char first_base, previous_base;  // Base at previous 2 positions
                if (query_pos >= 1) previous_base = seq_nt16_str[bam_seqi(read_seq, query_pos - 1)];
                else previous_base = 'N';
                if (query_pos >= 2) first_base = seq_nt16_str[bam_seqi(read_seq, query_pos - 2)];
                else first_base = 'N';
                char next_base, last_base;  // Base at next 2 positions
                if (query_pos + 1 < b->core.l_qseq) next_base = seq_nt16_str[bam_seqi(read_seq, query_pos + 1)];
                else next_base = 'N';
                if (query_pos + 2 < b->core.l_qseq) last_base = seq_nt16_str[bam_seqi(read_seq, query_pos + 2)];
                else last_base = 'N';

                // Check that bases indeed match or mismatch
                base_info.is_match = 1;

                // Check if CpG, CHH, CHG site
                Context context = ref.get_seq_context(first_base, previous_base, base, next_base, last_base, params->c_conversion, ref_pos, strand);
                base_info.context = context;

                if (params->c_conversion) {  // If unmethylated cytosines have been converted (either chemically or enzymatically)

                    // Check methylation status if CpG site, and check the presence of a putative mutation otherwise.
                    if ((is_cytosine(context)) && ((strand & 1) == 1)) {  // OT or CTOT case
                        assert(ref_base == 'C');
                        if (base == 'C') {
                            base_info.methylated = 1;
                        } else if (base == 'T') {
                            base_info.methylated = 0;
                        }
                    } else if ((is_cytosine(context)) && ((strand & 1) == 0)) {  // OB or CTOB
                        assert(ref_base == 'G');
                        if (base == 'G') {
                            base_info.methylated = 1;
                        } else if (base == 'A') {
                            base_info.methylated = 0;
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
                            } else {
                                // Ambiguous case: a mismatch may reflect bisulfite conversion OR a true
                                // mutation. In cfDNA and especially cancer-derived cfDNA, assuming every
                                // such event is a converted cytosine biases methylation calls toward
                                // unmethylated statuses. Keep the base conservative instead of forcing a
                                // methylation call in this scenario.
                                base_info.putative_snp = 1;
                                base_info.uninformative = 1;
                            }
                        } else {  // Base match
                            base_info.uninformative = 0;
                        }
                    }
                } else {  // No conversion of unmethylated cytosines
                    // In the absence of bisulfite/enzymatic conversion, any base mismatch
                    // is treated as a putative point mutation (or sequencing error)
                    if (base != ref_base) base_info.putative_snp = 1;
                }

                // If reference base is unknown, set it to the observed one
                if (ref_base_info.c == 'N') ref_base_info.c = base;

            } else if (cigar_op_type == BAM_CINS) {  // I (insertion to the reference)
                char base = seq_nt16_str[bam_seqi(read_seq, query_pos)];
                base_info.c = base;
                base_info.phred_score = read_qual[query_pos];
                base_info.inserted = 1;
                base_info.uninformative = 1;  // No reference -> no way to disentangle mutations from sequencing errors
                ref_base_info.deleted = true;
                consume_ref = false;
            } else if (cigar_op_type == BAM_CDEL) {  // D (deletion from the reference)
                base_info.uninformative = 1;
                base_info.deleted = 1;
                consume_query = false;
                valid_base = false;
            } else if (cigar_op_type == BAM_CREF_SKIP) {  // N (skipped region from the reference)
                base_info.uninformative = 1;
                consume_query = false;
                valid_base = false;
            } else if ((cigar_op_type == BAM_CSOFT_CLIP) || (cigar_op_type == BAM_CBACK)) {  // S (soft clipping) or B ("backwards" skip operator)
                // Soft-clipped base: base that is not part of the alignment
                base_info.uninformative = 1;
                base_info.soft_clipped = 1;
                char base = seq_nt16_str[bam_seqi(read_seq, query_pos)];
                base_info.c = base;
                base_info.phred_score = read_qual[query_pos];
                ref_base_info.soft_clipped = 1;
                consume_ref = false;
                valid_base = false;
            } else {  // Ignore H (hard clipping) and P (padding) operations
                // Hard-clipped base: base that is not part of the alignment and
                // is not present in the sequence at all.
                consume_query = false;
                consume_ref = false;
                valid_base = false;
            }

            // If base is unknown, then it is considered as non-informative
            if (base_info.c == 'N') base_info.uninformative = 1;

            // If Phred score is too low, mark base as uninformative
            if (base_info.phred_score < params->min_phred) base_info.uninformative = 1;

            // Store the alignment
            if (!(base_info.soft_clipped || ref_base_info.soft_clipped)) {
                if (base_info.inserted) read_info.add_base(ref_pos, base_info);
                else read_info.add_base(ref_pos + 1, base_info);
                read_info.add_ref_base(ref_pos + 1, ref_base_info);
            }

            // Go to next base
            if (consume_ref) ref_pos++;
            if (consume_query) query_pos++;
        }
    }

    // Detect unreliable bases at the ends of the read
    read_info.trim_alignment_inplace(*params);

    // Compute average Phred score
    uint32_t avg_phred = 0;
    for (size_t pos = 0; pos < b->core.l_qseq; pos++) {
        avg_phred += read_qual[pos];
    }
    avg_phred /= b->core.l_qseq;
    read_info.avg_phred = avg_phred;

    // Get read weight
    read_info.weight = 1.0f;
    uint8_t* gc_tag = bam_aux_get(b, "GC");
    if (gc_tag) {
        read_info.weight = bam_aux2f(gc_tag);
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

    kstring_t* c_os = (kstring_t*) calloc(1, sizeof(kstring_t));
    if (!c_os) {
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

    // Queue used to temporarily store a read's mate
    OrderedQueueMap<std::string, ReadInfo> paired_reads(10000);
    std::unordered_set<std::string> blacklisted_pairs;

    // Queue used to store methylation statuses
    OrderedQueueMap<genome_size_t, std::tuple<uint8_t, uint16_t, uint16_t, uint8_t>> meth_queue(100000);

    size_t n_processed_reads = 0;

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
            // If we've specified a region, then break once we're outside of it
            globalTid = std::numeric_limits<uint32_t>::max();
        }
        if (localTid < hdr->n_targets && (globalTid != std::numeric_limits<uint32_t>::max())) {
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
                id,
                hdr->target_name[b->core.tid],
                b->core.pos + 1,  // Change 0-based to 1-based coordinate
                bam_endpos(b),  // 0-based position of first base after the alignment
                abs(b->core.isize),  // Insert size
                get_strand(b),
                b->core.qual,
                b->core.flag
            );
            n_processed_reads++;

            // If secondary alignment, or mapping quality is insufficient, discard read.
            if ((!config->single_read) and (blacklisted_pairs.find(read.id) != blacklisted_pairs.end())) {
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

            // If read is not part of allowed read groups, skip it.
            if (config->use_read_group) {
                uint8_t* rg_tag = bam_aux_get(b, "RG");
                if (rg_tag) {
                    const char* rg_value = bam_aux2Z(rg_tag);
                    if (!std::regex_match(rg_value, config->read_group)) continue;
                }
            }

            // Use inclusion bounds to trim alignments
            trim_alignment(b, *config);

            // Extract read-level statistics from read
            AlignmentRefSequence ref(seq, localPos2, seq_len);
            process_read(b, ref, read, config);

            // Call methylation using MM and ML tags (if required)
            if (config->use_mm_ml) methylation_from_aux(b, read, config->ml_cutoff);
            // TODO: call methylation from aligner-specific optional fields such as Bismark (XM tag), if using them

            // Combine read pairs into fragments
            if (config->single_read) {  // Single-read
                Fragment fragment = compute_fragment_statistics(read);
                if (fragment.is_valid() && (fragment.length <= config->max_fragment_size)) {
                    if (config->write_msi) detect_microsatellites_in_read(config->msi_output_fp, read, config->msi_flank, fragment.weight);
                    write_fragment(os, fragment);
                    if (config->write_cytosines) write_methylation(c_os, meth_queue, fragment);
                }
            } else {  // Paired-end
                bool found_pair = false;
                if (paired_reads.contains(read.id)) {
                    ReadInfo& other_read = paired_reads.get(read.id);
                    if (read.paired_with(other_read)) {
                        Fragment fragment = compute_fragment_statistics(read, other_read);
                        if (fragment.is_valid() && (read.chr_id == other_read.chr_id) && (fragment.length <= config->max_fragment_size)) {
                            if (config->write_msi) detect_microsatellites_in_pair(config->msi_output_fp, read, other_read, config->msi_flank, fragment.weight);
                            write_fragment(os, fragment);
                            if (config->write_cytosines) write_methylation(c_os, meth_queue, fragment);
                        }
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
            if (config->write_cytosines) {
                if (c_os->l) fwrite(c_os->s, 1, c_os->l, config->c_output_fp[0]);
                c_os->l = 0;
            }
            outputBin++;
            pthread_mutex_unlock(&outputMutex);
            break;
        }
    }

    if (os->s) free(os->s);
    free(os);
    if (c_os->s) free(c_os->s);
    free(c_os);
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
" -m STR     Write microsatellite calls to the given text file (4 columns).\n"
" -x INT     Minimum number of flanking bases required around a microsatellite\n"
"            for MSI detection (default 3).\n"
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
    FILE *c_ofile = stdout;
    faidx_t *fai;
    std::string out_filepath;

    static struct option lopts[] = {
            {"help",    0, nullptr, 'h'},
            {"version", 0, nullptr, 'v'},
            {"single_read",   0, nullptr, 20},
            {"useMMML",   0, nullptr, 21},
            {"mlCutoff",   1, nullptr, 22},
            {"cconv",   0, nullptr, 23},
            {"max_frag_size",   1, nullptr, 24},
            {"read_group",   1, nullptr, 25},
            {"partition_size",    1, nullptr,  19},
            {"msi-output",   1, nullptr, 26},
            {"msi-flank",   1, nullptr, 27},
            {"OT",           1, nullptr,   7},
            {"OB",           1, nullptr,   8},
            {"CTOT",         1, nullptr,   9},
            {"CTOB",         1, nullptr,  10},
            {0,         0, nullptr,   0}
    };
    char c;
    while((c = getopt_long(argc, argv, "hvq:p:o:m:n:c:@:r:l:F:R:x:", lopts, nullptr)) >= 0) {
        switch(c) {
            case 'h':
                cli_help();
                return 0;
            case 'o':
                out_filepath = optarg;
                if ((ofile = fopen(optarg, "wb")) == nullptr) {
                    fprintf(stderr, "Couldn't open %s for writing\n", optarg);
                    return 2;
                }
                break;
            case 'm':
                config.write_msi = true;
                if ((config.msi_output_fp = fopen(optarg, "w")) == nullptr) {
                    fprintf(stderr, "Couldn't open %s for writing\n", optarg);
                    return 2;
                }
                break;
            case 'x':
                config.msi_flank = static_cast<uint32_t>(std::max(0, atoi(optarg)));
                break;
            case 'c':
                if ((c_ofile = fopen(optarg, "wb")) == nullptr) {
                    fprintf(stderr, "Couldn't open %s for writing\n", optarg);
                    return 2;
                }
                config.write_cytosines = true;
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
            case 23:
                config.c_conversion = 1;
                break;
            case 25:
                config.read_group = std::regex(optarg);
                config.use_read_group = true;
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
            case 24:
                config.max_fragment_size = strtoul(optarg, nullptr, 10);
                break;
            case 19:
                config.partition_size = strtoul(optarg, nullptr, 10);
                if(config.partition_size < 1) {
                    fprintf(stderr, "Error: The chunk size must be at least 1!\n");
                    return 1;
                }
                break;
            case 26:
                config.write_msi = true;
                if ((config.msi_output_fp = fopen(optarg, "w")) == nullptr) {
                    fprintf(stderr, "Couldn't open %s for writing\n", optarg);
                    return 2;
                }
                break;
            case 27:
                config.msi_flank = static_cast<uint32_t>(std::max(0, atoi(optarg)));
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
    if (config.msi_flank > 10000) {
        fprintf(stderr, "-x %u is unusually large; clamping it to 10000.\n", config.msi_flank);
        config.msi_flank = 10000;
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

    // MSI output has chromosome ID so coordinates are not ambiguous.
    if (config.write_msi && config.msi_output_fp) {
        fprintf(
            config.msi_output_fp,
            "#chrom_id\tstart0\tend0\tobserved_length\tmotif\t"
            "reference_length\tobserved_repeat_count\t"
            "reference_repeat_count\tfragment_weight\n"
        );
    }

    // Output files
    config.output_fp = (FILE**) malloc(sizeof(FILE *));
    config.output_fp[0] = ofile;
    config.c_output_fp = (FILE**) malloc(sizeof(FILE *));
    config.c_output_fp[0] = c_ofile;

    // Spin up the threads
    pthread_mutex_init(&positionMutex, nullptr);
    pthread_mutex_init(&outputMutex, nullptr);
    pthread_t *threads = (pthread_t*) calloc(config.n_threads, sizeof(pthread_t));
    for (int i = 0; i < config.n_threads; i++) pthread_create(threads + i, nullptr, &process_partition, &config);
    for (int i = 0; i < config.n_threads; i++) pthread_join(threads[i], nullptr);
    free(threads);
    pthread_mutex_destroy(&positionMutex);
    pthread_mutex_destroy(&outputMutex);
    if (config.write_msi) {
        fprintf(stderr, "MSI: reads searched=%" PRIu64
                ", candidate repeats=%" PRIu64 ", records written=%" PRIu64 "\n",
                msi_reads_checked.load(), msi_candidates_found.load(),
                msi_records_written.load());
    }

    // Add dummy fragment to output file to indicate that the file is valid
    Fragment dummy;
    dummy.chr_id = 0;
    dummy.length = 0;
    char buffer[FRAGMENT_SIZE_IN_BYTES];
    dummy.serialize(buffer);
    fwrite(buffer, sizeof(buffer), 1, config.output_fp[0]);

    // Save overall methylation counts in separate file
    write_methylation_counts_to_file(out_filepath + ".jag");

    // Close files and free memory
    hts_close(config.fp);
    if (ofile != stdout) fclose(ofile);
    if (c_ofile != stdout) fclose(c_ofile);
    if (config.msi_output_fp != nullptr) fclose(config.msi_output_fp);
    free(config.output_fp);
    free(config.c_output_fp);
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
