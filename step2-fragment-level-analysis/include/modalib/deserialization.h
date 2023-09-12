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

#endif // MODALIB_DESERIALIZATION_H
