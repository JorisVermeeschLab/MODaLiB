#ifndef _SERIALIZATION_H__
#define _SERIALIZATION_H__

#include <iostream>
#include <fstream>
#include <string>


class Serializer {
public:
    inline static int serialize(char* buffer, std::uint8_t value) {
        return Serializer::serialize_int(buffer, value);
    }
    inline static int serialize(char* buffer, std::uint16_t value) {
        return Serializer::serialize_int(buffer, value);
    }
    inline static int serialize(char* buffer, std::uint32_t value) {
        return Serializer::serialize_int(buffer, value);
    }
private:
    template <typename Int>
    static int serialize_int(char* buffer, Int value) {
        typedef typename std::make_unsigned<Int>::type UInt;
        constexpr UInt byte{0xFF};
        for(std::size_t i{0}; i < sizeof(Int); ++i)
            buffer[i] = static_cast<char>(
                ((static_cast<UInt>(value) & (byte << (8*i))) >> (8*i)) & 0xFF
            );
        return sizeof(Int);
    }
};

#endif // _SERIALIZATION_H__
