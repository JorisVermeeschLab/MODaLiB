#ifndef SERIALIZATION_H__
#define SERIALIZATION_H__

#include <iostream>
#include <fstream>
#include <string>


class Serializer {
public:
    inline static int serialize(char* buffer, uint8_t value) {
        return Serializer::serialize_int(buffer, value);
    }
    inline static int serialize(char* buffer, uint16_t value) {
        return Serializer::serialize_int(buffer, value);
    }
    inline static int serialize(char* buffer, uint32_t value) {
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

#endif // SERIALIZATION_H__
