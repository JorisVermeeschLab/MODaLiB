//  sample.cpp
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

#include "modalib/sample.h"
#include "modalib/utils.h"
#include "modalib/deserialization.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


namespace modalib {

_Sample::_Sample() {
	// Nothing to do, as the location of the input file is unknown
}

_Sample::_Sample(const _Sample& other) : 
        cursor(other.cursor - 1),
        lb(other.lb),
        ub(other.ub),
        filepath(other.filepath),
        frag(other.frag),
        previous_frag(other.previous_frag) {
    this->open(other.filepath);
    this->go_to_fragment(other.cursor - 1);
}

_Sample::_Sample(std::string &filepath) : 
        cursor(0),
        lb(0),
        ub(std::numeric_limits<long>::max()),
        filepath(filepath),
        frag(),
        previous_frag() {
    this->open(filepath);
    this->go_to_fragment(0);
}

void _Sample::operator=(const _Sample& other) {
    this->cursor = other.cursor - 1;
    this->lb = other.lb;
    this->ub = other.ub;
    this->filepath = other.filepath;
    this->open(other.filepath);
    this->go_to_fragment(other.cursor - 1);
}

_Sample::~_Sample() {
	this->is.close();
}

void _Sample::open(const std::string &filepath) {
    if (!this->is.is_open()) {
        this->is.open(filepath, std::ios::in | std::ios::binary);
    }
    if (!this->is.is_open()) throw std::runtime_error("Could not open file " + filepath);
    this->is.seekg(0, std::iostream::beg);
}

void _Sample::next() {

    // Check if end-of-file has been reached
	if (this->is.eof()) return;

	// Keep track of previous frag
	this->previous_frag = this->frag;

    // Load new frag
    char buffer[FRAGMENT_SIZE_IN_BYTES];
    this->is.read(buffer, FRAGMENT_SIZE_IN_BYTES);

    this->cursor++;

    this->frag.deserialize(buffer);
}

void _Sample::go_to_fragment(long i) {
    if (i <= 0) {
        this->cursor = 0;
        this->is.seekg(0, std::iostream::beg);
        this->next();
    } else {
        this->cursor = i - 1;
        this->is.seekg((i - 1) * FRAGMENT_SIZE_IN_BYTES, std::iostream::beg);
        this->next();
        this->next();
    }
}

bool _Sample::eof() {
    return ((this->cursor >= this->ub) || (this->is.eof()));
}

void _Sample::set_lower_bound(long n) {
    this->lb = n;
    this->go_to_fragment(n);
}

void _Sample::set_upper_bound(long n) {
    this->ub = n;
}

} // namespace modalib
