//  sample.h
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

#ifndef MODALIB_SAMPLE_H
#define MODALIB_SAMPLE_H

#include "fragment_info.h"
#include "methylation.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace modalib {

class _Sample {
private:

    // Index of the current fragment
	long cursor = 0;

    // Interval of the fragments to iterate over
    long lb = 0;
	long ub = std::numeric_limits<long>::max();

    // Path to the input binary file
	std::string filepath;

    // Input stream
	std::ifstream is;

	// Whether input file contains fragments or cytosine contexts
	bool contains_fragments = true;

    // Last fragment parsed
	FragmentInfo frag;

    // Second-last fragment parsed
	FragmentInfo previous_frag;

	// Last cytosine context parsed
	CytosineContext c;

    // Second-last cytosine context parsed
    CytosineContext previous_c;

public:
	_Sample() = default;
	_Sample(const _Sample&);
	explicit _Sample(const std::string &filepath);
	_Sample& operator=(const _Sample&);
    ~_Sample();
    void open(const std::string &filepath);
    void go_to_chunk(long i);
	void next();
	bool eof() const;

    // Setters
    void set_lower_bound(long n);
    void set_upper_bound(long n);

    // Getters and setters
    inline const FragmentInfo& get() const {
        if (!contains_fragments) throw std::runtime_error("Method unavailable for .cytosines.modalib files");
        return frag;
    }
    inline const FragmentInfo& get_previous() const {
        if (!contains_fragments) throw std::runtime_error("Method unavailable for .cytosines.modalib files");
        return previous_frag;
    }
    inline const CytosineContext& get_c() const {
        if (contains_fragments) throw std::runtime_error("Method unavailable for .modalib files");
        return c;
    }
    inline const CytosineContext& get_previous_c() const {
        if (contains_fragments) throw std::runtime_error("Method unavailable for .modalib files");
        return previous_c;
    }
    inline long get_cursor() const { return cursor; }

    long get_length();
private:
    size_t get_chunk_size() const { return (contains_fragments) ? FRAGMENT_SIZE_IN_BYTES : CYTOSINE_CONTEXT_SIZE_IN_BYTES; }
};

} // namespace modalib

#endif  // MODALIB_SAMPLE_H
