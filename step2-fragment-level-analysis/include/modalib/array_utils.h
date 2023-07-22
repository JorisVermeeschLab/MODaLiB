#ifndef _ARRAY_UTILS_H__
#define _ARRAY_UTILS_H__

#include <cmath>
#include <valarray>

namespace modalib {

inline std::valarray<long> array_max(const std::valarray<long>& arr, long value) {
    std::valarray<long> res(arr.size());
    for (size_t i = 0; i < arr.size(); i++) {
        res[i] = std::max(arr[i], value);
    }
    return res;
}

inline int argmin(const std::valarray<long>& coords) noexcept {
    const int n = coords.size();
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (coords[i] < coords[k]) k = i;
    }
    return k;
}

} // namespace modalib

#endif // _ARRAY_UTILS_H__
