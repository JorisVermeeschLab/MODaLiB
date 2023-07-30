#ifndef _SITE_H__
#define _SITE_H__

#include <cmath>
#include <functional>
#include <stdint.h>

#include "fragment_info.h"

namespace modalib {

struct Site {
    uint8_t chr;
    uint32_t pos;
};

class Region {
public:
    Site start;
    Site end;
    Region(uint8_t chr, uint32_t start, uint32_t end);
    int cmp(const FragmentInfo&);
    bool operator==(const Region &other) const {
        bool v = (start.chr == other.start.chr);
        v &= (start.pos == other.start.pos);
        v &= (end.chr == other.end.chr);
        v &= (end.pos == other.end.pos);
        return v;
    }
};

struct RegionHashFunction {
    std::hash<uint8_t> hash_func8;
    std::hash<uint32_t> hash_func32;
    // TODO: optimize hash function: CRC32 hash algorithm?
    std::size_t operator()(const Region &r) const {
        std::size_t h = 157;  // TODO
        h = h ^ hash_func8(r.start.chr);
        h = h ^ hash_func32(r.start.pos);
        h = h ^ hash_func8(r.end.chr);
        h = h ^ hash_func32(r.end.pos);
        return h;
    }
};

} // namespace modalib

#endif // _SITE_H__