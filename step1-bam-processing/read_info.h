#ifndef _READ_INFO_H__
#define _READ_INFO_H__

#include <string.h>
#include <cstdlib>


inline uint8_t parse_chr_name(const char* c) {
    size_t n = strlen(c);
    uint8_t chr;
    if ((n > 3) && (c[0] == 'c')) {
        c = c + 3;
    }
    char *p;
    long converted = strtol(c, &p, 10);
    if (*p) {
        if (n > 3) {
            if (c[0] == 'X') chr = 23;
            else if (c[0] == 'Y') chr = 24;
            else if (c[0] == 'M') chr = 25;
            else chr = 0;
        } else if ((n == 3) && (strncmp(c, "EBV", 3) == 0)) {
            chr = 26;
        } else {
            chr = 0;
        }
    }
    else {
        chr = (uint8_t) converted;
    }
    return chr;
}


struct ReadInfo {
public:

    // Read identifier
    std::string id;

    // Chromosome id
    uint8_t chr_id;

    // Left-most coordinate on chromosome
    uint32_t start;

    // Number of methylated CpG sites
    uint16_t n_methylated = 0;
    uint16_t n_methylated_nov = 0;  // Non-overlapping

    // Number of non-methylated CpG sites
    uint16_t n_unmethylated = 0;
    uint16_t n_unmethylated_nov = 0;  // Non-overlapping

    // Number of single nucleotide polymorphisms
    uint16_t snp_count = 0;
    uint16_t snp_count_nov = 0;  // Non-overlapping

    // Number of cytosines and guanines
    uint16_t gc_count = 0;
    uint16_t gc_count_nov = 0;  // Non-overlapping

    // Number of informative bases
    uint16_t s_count = 0;
    uint16_t s_count_nov = 0;  // Non-overlapping

    // Length of the query sequence (read)
    uint16_t length;

    // Insert size
    uint16_t insert_size;

    // SAM flag
    uint16_t sam_flag;

    // Number of inserted bases
    uint16_t inserted = 0;

    // Number of deleted bases
    uint16_t deleted = 0;

public:
    ReadInfo() = default;
    ReadInfo(std::string& id, const char* chr_name, uint32_t start, uint16_t length, uint16_t insert_size, uint16_t sam_flag) :
        id(id),
        chr_id(parse_chr_name(chr_name)),
        start(start),
        length(length),
        insert_size(insert_size),
        sam_flag(sam_flag) {}
    ReadInfo(const ReadInfo& other) = default;
    bool is_first_in_pair() { return 0x40 & sam_flag; }
    bool is_reverse_strand() { return 0x10 & sam_flag; }
    bool mapped_in_proper_pair() { return 0x3 & sam_flag; }
    bool failed_quality_checks() { return 0x200 & sam_flag; }
    bool is_duplicate() { return 0x400 & sam_flag; }
    bool paired_with(ReadInfo& other) { return this->id == other.id; }
    bool is_non_overlapping(uint32_t read_position) {
        if (!is_first_in_pair()) return true;
        else {
            // TODO: We assume that the 2 paired reads have the same length
            // TODO: this simplifies the computation of the overlap
            return (int) read_position < (int) insert_size - (int) length;
        }
    }
    void add_base(int base, uint32_t read_position, int strand) {
        this->add_base(base, read_position, strand, false);
    }
    void add_base(int base, uint32_t read_position, int strand, bool is_snp) {
        if (is_snp) inc_snp_count(read_position);
        inc_s_count(read_position);
        if ((strand & 1) == 1) {
            if (base == 2) inc_gc_count(read_position);  // C
            else if (base == 4) inc_gc_count(read_position);  // G
        } else if ((strand & 1) == 0) {
            if (base == 1) inc_gc_count(read_position);  // A
            else if (base == 8) inc_gc_count(read_position);  // T
        }
    }
    void inc_gc_count(uint32_t read_position) {
        gc_count++;
        if (is_non_overlapping(read_position)) gc_count_nov++;
    }
    void inc_s_count(uint32_t read_position) {
        s_count++;
        if (is_non_overlapping(read_position)) s_count_nov++;
    }
    void inc_methyl_count(uint32_t read_position) {
        n_methylated++;
        if (is_non_overlapping(read_position)) n_methylated_nov++;
    }
    void inc_unmethyl_count(uint32_t read_position) {
        n_unmethylated++;
        if (is_non_overlapping(read_position)) n_unmethylated_nov++;
    }
    void inc_snp_count(uint32_t read_position) {
        snp_count++;
        if (is_non_overlapping(read_position)) snp_count_nov++;
    }
};

#endif // _READ_INFO_H__
