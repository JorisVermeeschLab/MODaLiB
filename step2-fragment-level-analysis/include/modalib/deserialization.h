//  deserialization.h
//
//  Copyright 2023 Antoine Passemiers <antoine.passemiers@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//  MA 02110-1301, USA.

#ifndef MODALIB_DESERIALIZATION_H
#define MODALIB_DESERIALIZATION_H

#include <cstring>
#include <string>
#include <vector>

namespace modalib {


class Deserializer {
public:
    inline static float int_as_float(std::uint32_t input) {
        float output;
        char* ip = (char*) &input, *op = (char*) &output;
        std::memcpy(op, ip, sizeof(input));
        return output;
    }
    inline static int serialize(char* buffer, const std::uint8_t& input) {
        return Deserializer::serialize_int(buffer, input);
    }
    inline static int serialize(char* buffer, const std::uint16_t& input) {
        return Deserializer::serialize_int(buffer, input);
    }
    inline static int serialize(char* buffer, const std::uint32_t& input) {
        return Deserializer::serialize_int(buffer, input);
    }
    inline static int serialize(char* buffer, float input) {
        uint32_t ivalue = *((uint32_t*) &input);
        return Deserializer::serialize_int(buffer, ivalue);
    }
    inline static int deserialize(const char* buffer, std::uint8_t& output) {
        return Deserializer::deserialize_int(buffer, output);
    }
    inline static int deserialize(const char* buffer, std::uint16_t& output) {
        return Deserializer::deserialize_int(buffer, output);
    }
    inline static int deserialize(const char* buffer, std::uint32_t& output) {
        return Deserializer::deserialize_int(buffer, output);
    }
    inline static int deserialize(const char* buffer, float& output) {
        uint32_t temp;
        int offset = Deserializer::deserialize_int(buffer, temp);
        output = int_as_float(temp);
        return offset;
    }
private:
    template <typename Int>
    static int serialize_int(char* buffer, const Int& input) {
        Int v = input;
        typedef typename std::make_unsigned<Int>::type UInt;
        for (std::size_t i{0}; i < sizeof(Int); ++i) {
            buffer[i] = static_cast<char>((input >> (i * 8)) & 0xFF);
        }
        return sizeof(Int);
    }
    template <typename Int>
    static int deserialize_int(const char* buffer, Int& output) {
        typedef typename std::make_unsigned<Int>::type UInt;
        // constexpr UInt byte{0xFF};
        output = 0;
        for (std::size_t i{0}; i < sizeof(Int); ++i) {
            output |= static_cast<UInt>(static_cast<unsigned char>(buffer[i])) << (i * 8);
        }
        return sizeof(Int);
    }
};


} // namespace modalib

#endif // MODALIB_DESERIALIZATION_H
