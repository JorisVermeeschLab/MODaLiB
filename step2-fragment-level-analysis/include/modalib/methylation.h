//  methylation.h
//
//  Copyright 2024 Antoine Passemiers <antoine.passemiers@gmail.com>
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

#ifndef MODALIB_METHYLATION_H
#define MODALIB_METHYLATION_H

#include <cassert>
#include <cstdint>
#include "deserialization.h"

namespace modalib {

// Number of bytes used to represent each cytosine context in the input binary file
static constexpr size_t CYTOSINE_CONTEXT_SIZE_IN_BYTES = 11;

struct CytosineContext {
public:
    uint8_t chr_id;
    uint32_t pos;
    uint8_t type;
    uint16_t n_methylated;
    uint16_t n_unmethylated;
    uint8_t strand;
    float meth_ratio() const {
        return static_cast<float>(this->n_methylated) / static_cast<float>(this->n_methylated + this->n_unmethylated);
    }
    void deserialize(const char* buffer) {
        int offset = 0;
        offset += Deserializer::deserialize(buffer + offset, this->chr_id);
        offset += Deserializer::deserialize(buffer + offset, this->pos);
        offset += Deserializer::deserialize(buffer + offset, this->type);
        offset += Deserializer::deserialize(buffer + offset, this->n_methylated);
        offset += Deserializer::deserialize(buffer + offset, this->n_unmethylated);
        offset += Deserializer::deserialize(buffer + offset, this->strand);
        assert(offset == CYTOSINE_CONTEXT_SIZE_IN_BYTES);
    }
};

}  // namespace modalib

#endif  // MODALIB_METHYLATION_H
