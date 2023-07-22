#ifndef _DESERIALIZATION_H__
#define _DESERIALIZATION_H__

#include <string>
#include <vector>

namespace modalib {


class Deserializer {
public:
    inline static int deserialize(const char* buffer, std::uint8_t& output) {
        return Deserializer::deserialize_int(buffer, output);
    }
    inline static int deserialize(const char* buffer, std::uint16_t& output) {
        return Deserializer::deserialize_int(buffer, output);
    }
    inline static int deserialize(const char* buffer, std::uint32_t& output) {
        return Deserializer::deserialize_int(buffer, output);
    }
private:
    template <typename Int>
    static int deserialize_int(const char* buffer, Int& output) {
        typedef typename std::make_unsigned<Int>::type UInt;
        // constexpr UInt byte{0xFF};
        output = 0;
        for (std::size_t i{0}; i < sizeof(Int); ++i) {
            output |= static_cast<UInt>(static_cast<unsigned char>(buffer[i])) << (8*i);
        }
        return sizeof(Int);
    }
};


} // namespace modalib

#endif // _DESERIALIZATION_H__
