//  sample.h
//
//  Copyright 2022 Antoine Passemiers <antoine.passemiers@gmail.com>
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

#ifndef _SAMPLE_H__
#define _SAMPLE_H__

#include "read_info.h"
#include "osc.h"

#include <fstream>
#include <iostream>
#include <limits>
#include <string>

namespace modalib {

// Constant indicating that EOF has been reached
const long END_OF_SAMPLE = std::numeric_limits<long>::max();

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

    // Last fragment parsed
	ReadInfo read_info;

    // Second-last fragment parsed
	ReadInfo previous_read;

public:
	_Sample();
	_Sample(const _Sample&);
	explicit _Sample(std::string &filepath);
	void operator=(const _Sample&);
    ~_Sample();
    void open(const std::string &filepath);
    void reset();
    void go_to_fragment(int i);
	void next();
	bool eof();

    // Setters
    void set_lower_bound(long n);
    void set_upper_bound(long n);

    // Getters and setters
    inline const ReadInfo& get() { return read_info; }
    inline const ReadInfo& get_previous() { return previous_read; }

};

} // namespace modalib

#endif // _SAMPLE_H__
