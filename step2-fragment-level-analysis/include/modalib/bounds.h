#ifndef _BOUNDS_H
#define _BOUNDS_H

#include <algorithm>
#include <string>
#include <iostream>
#include <vector>

namespace modalib {

static const long NUM_50KB_BINS = 61625L;
static const long BIN_SIZE = 50000L;
static const long NUM_PARTITIONS = 10L;

static const std::vector<std::pair<std::string, uint32_t>> CHROMOSOME_BOUNDS {
    std::pair<std::string, uint32_t>("Unknown", 0),
    std::pair<std::string, uint32_t>("chr1", 4945),
    std::pair<std::string, uint32_t>("chr2", 4860),
    std::pair<std::string, uint32_t>("chr3", 3991),
    std::pair<std::string, uint32_t>("chr4", 3826),
    std::pair<std::string, uint32_t>("chr5", 3618),
    std::pair<std::string, uint32_t>("chr6", 3418),
    std::pair<std::string, uint32_t>("chr7", 3177),
    std::pair<std::string, uint32_t>("chr8", 2926),
    std::pair<std::string, uint32_t>("chr9", 2806),
    std::pair<std::string, uint32_t>("chr10", 2708),
    std::pair<std::string, uint32_t>("chr11", 2690),
    std::pair<std::string, uint32_t>("chr12", 2647),
    std::pair<std::string, uint32_t>("chr13", 2283),
    std::pair<std::string, uint32_t>("chr14", 2128),
    std::pair<std::string, uint32_t>("chr15", 2007),
    std::pair<std::string, uint32_t>("chr16", 1777),
    std::pair<std::string, uint32_t>("chr17", 1576),
    std::pair<std::string, uint32_t>("chr18", 1523),
    std::pair<std::string, uint32_t>("chr19", 1277),
    std::pair<std::string, uint32_t>("chr20", 1249),
    std::pair<std::string, uint32_t>("chr21", 939),
    std::pair<std::string, uint32_t>("chr22", 994),
    std::pair<std::string, uint32_t>("chrX", 3099),
    std::pair<std::string, uint32_t>("chrY", 1156),
    std::pair<std::string, uint32_t>("chrM", 1),
    std::pair<std::string, uint32_t>("EBV", 4)
};

static const std::vector<std::pair<long, long>> CHROMOSOME_BOUNDS_CUMULATIVE {
    std::pair<long, long>(0L, 0L),
    std::pair<long, long>(0L, 4945L),
    std::pair<long, long>(4945L, 4860L),
    std::pair<long, long>(9805L, 3991L),
    std::pair<long, long>(13796L, 3826L),
    std::pair<long, long>(17622L, 3618L),
    std::pair<long, long>(21240L, 3418L),
    std::pair<long, long>(24658L, 3177L),
    std::pair<long, long>(27835L, 2926L),
    std::pair<long, long>(30761L, 2806L),
    std::pair<long, long>(33567L, 2708L),
    std::pair<long, long>(36275L, 2690L),
    std::pair<long, long>(38965L, 2647L),
    std::pair<long, long>(41612L, 2283L),
    std::pair<long, long>(43895L, 2128L),
    std::pair<long, long>(46023L, 2007L),
    std::pair<long, long>(48030L, 1777L),
    std::pair<long, long>(49807L, 1576L),
    std::pair<long, long>(51383L, 1523L),
    std::pair<long, long>(52906L, 1277L),
    std::pair<long, long>(54183L, 1249L),
    std::pair<long, long>(55432L, 939L),
    std::pair<long, long>(56371L, 994L),
    std::pair<long, long>(57365L, 3099L),
    std::pair<long, long>(60464L, 1156L),
    std::pair<long, long>(61620L, 1L),
    std::pair<long, long>(61621L, 4L)
};

inline unsigned long clamp(unsigned long value, unsigned long lb, unsigned long ub) {
    if (value < lb) value = lb;
    if (value > ub) value = ub;
    return value;
}

inline int coords_to_partition(uint8_t chr_id, long pos, size_t n_partitions) {
    const std::pair<long, long>& coords = CHROMOSOME_BOUNDS_CUMULATIVE.at(chr_id);
    const unsigned long combined_coords = coords.first * BIN_SIZE + clamp(pos, 0L, coords.second * BIN_SIZE - 1L);
    size_t partition_id = static_cast<unsigned long>(combined_coords * n_partitions) / static_cast<unsigned long>(BIN_SIZE * NUM_50KB_BINS);
    return clamp(partition_id, 0, n_partitions - 1);
}

} // namespace modalib

#endif // _BOUNDS_H
