#ifndef _BOUNDS_H
#define _BOUNDS_H

#include <algorithm>
#include <string>
#include <vector>

namespace modalib {

static const uint32_t NUM_50KB_BINS = 61625;
static const uint32_t BIN_SIZE = 50000;
static const uint32_t NUM_PARTITIONS = 10;

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

static const std::vector<std::pair<uint32_t, uint32_t>> CHROMOSOME_BOUNDS_CUMULATIVE {
    std::pair<uint32_t, uint32_t>(0, 0),
    std::pair<uint32_t, uint32_t>(0, 4945),
    std::pair<uint32_t, uint32_t>(4945, 4860),
    std::pair<uint32_t, uint32_t>(9805, 3991),
    std::pair<uint32_t, uint32_t>(13796, 3826),
    std::pair<uint32_t, uint32_t>(17622, 3618),
    std::pair<uint32_t, uint32_t>(21240, 3418),
    std::pair<uint32_t, uint32_t>(24658, 3177),
    std::pair<uint32_t, uint32_t>(27835, 2926),
    std::pair<uint32_t, uint32_t>(30761, 2806),
    std::pair<uint32_t, uint32_t>(33567, 2708),
    std::pair<uint32_t, uint32_t>(36275, 2690),
    std::pair<uint32_t, uint32_t>(38965, 2647),
    std::pair<uint32_t, uint32_t>(41612, 2283),
    std::pair<uint32_t, uint32_t>(43895, 2128),
    std::pair<uint32_t, uint32_t>(46023, 2007),
    std::pair<uint32_t, uint32_t>(48030, 1777),
    std::pair<uint32_t, uint32_t>(49807, 1576),
    std::pair<uint32_t, uint32_t>(51383, 1523),
    std::pair<uint32_t, uint32_t>(52906, 1277),
    std::pair<uint32_t, uint32_t>(54183, 1249),
    std::pair<uint32_t, uint32_t>(55432, 939),
    std::pair<uint32_t, uint32_t>(56371, 994),
    std::pair<uint32_t, uint32_t>(57365, 3099),
    std::pair<uint32_t, uint32_t>(60464, 1156),
    std::pair<uint32_t, uint32_t>(61620, 1),
    std::pair<uint32_t, uint32_t>(61621, 4)
};

inline int coords_to_partition(uint8_t chr_id, long pos, size_t n_partitions) {
    const std::pair<uint32_t, uint32_t>& coords = CHROMOSOME_BOUNDS_CUMULATIVE.at(chr_id);
    const long combined_coords = coords.first * BIN_SIZE + std::min(pos, static_cast<long>(coords.second * BIN_SIZE - 1));
    size_t partition_id = (combined_coords * n_partitions) / static_cast<long>(BIN_SIZE * NUM_50KB_BINS);
    return std::max(0, static_cast<int>(std::min(partition_id, n_partitions - 1)));
}

} // namespace modalib

#endif // _BOUNDS_H
