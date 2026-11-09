//  translocation.h
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

#ifndef MODALIB_TRANSLOCATION_H
#define MODALIB_TRANSLOCATION_H

#include <cstdint>
#include <utility>
#include <sstream>
#include <string>
#include <cmath>

#include "deserialization.h"
#include "sample.h"

namespace modalib {

void detect_putative_translocations(const _Sample& _sample) {
    _Sample sample(_sample);
    while (!sample.eof()) {
        const FragmentInfo& frag = sample.get();

        if (frag.chromosomes_mismatch()) {
            std::cout << 1*frag.start_chr_id << ", " << 1*frag.end_chr_id << std::endl;
        }

        sample.next();
    }
}

} // namespace modalib

#endif  // MODALIB_TRANSLOCATION_H
