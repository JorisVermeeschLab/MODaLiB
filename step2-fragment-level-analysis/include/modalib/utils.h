#ifndef _UTILS_H__
#define _UTILS_H__

#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

namespace modalib {

// Constant indicating that EOF has been reached
const long END_OF_SAMPLE = std::numeric_limits<long>::max();

typedef unsigned int count_t;

inline long combine_coords(uint8_t chr_id, long pos) {
    if (pos == END_OF_SAMPLE) return END_OF_SAMPLE;
    else return static_cast<long>(chr_id) * 1e+9 + pos;
}

inline int8_t parse_chr_name(const std::string &chr_name) {
    std::size_t n = chr_name.size();
    uint8_t chr;
    const char *c = chr_name.c_str();
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
        chr = static_cast<uint8_t>(converted);
    }
    return chr;
}

} // namespace modalib

#endif // _UTILS_H__
