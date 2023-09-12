//  reference.h
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

#ifndef MODALIB_REFERENCE_H
#define MODALIB_REFERENCE_H

#include <cctype>

struct AlignmentRefSequence {
private:
    char* seq;
    int start_pos;
    int length;
public:
    AlignmentRefSequence(char* seq, int start_pos, int length) {
        this->seq = seq;
        this->start_pos = start_pos;
        this->length = length;
    }
    ~AlignmentRefSequence() = default;  // Don't delete `seq` pointer
    char get_base(int abs_pos) const {
        return *(this->seq + abs_pos - this->start_pos);
    }
    std::string get_substring(int start, int end) const {
        return std::string(seq + start, end - start);
    }
    int is_cpg(int abs_pos) const {
        int pos = abs_pos - this->start_pos;
        if (pos >= this->length) return 0;
        char base = static_cast<char>(std::toupper(*(this->seq + pos)));
        if (base == 'C') {
            if (pos + 1 == this->length) return 0;
            char next_base = static_cast<char>(std::toupper(*(this->seq + pos + 1)));
            if (next_base == 'G') return 1;
            return 0;
        } else if (base == 'G') {
            if (pos == 0) return 0;
            char previous_base = static_cast<char>(std::toupper(*(this->seq + pos - 1)));
            if (previous_base == 'C') return -1;
            return 0;
        }
        return 0;
    }
    int is_ch(int ref_pos) const {
        int pos = ref_pos - this->start_pos;
        if (pos >= this->length) return 0;
        char base = static_cast<char>(std::toupper(*(this->seq + pos)));
        if (base == 'C') {
            if (pos + 1 == this->length) return 0;
            char next_base = static_cast<char>(std::toupper(*(this->seq + pos + 1)));
            if (next_base != 'G') return 1;
            return 0;
        } else if (base == 'G') {
            if (pos == 0) return 0;
            char previous_base = static_cast<char>(std::toupper(*(this->seq + pos - 1)));
            if (previous_base != 'C') return -1;
            return 0;
        }
        return 0;
    }
    int is_chg(int ref_pos) const {
        int pos = ref_pos - this->start_pos;
        if (pos >= this->length) return 0;
        char base = static_cast<char>(std::toupper(*(this->seq + pos)));
        if (base == 'C') {
            if (pos + 2 >= this->length) return 0;
            char second_base = static_cast<char>(std::toupper(*(this->seq + pos + 1)));
            char third_base = static_cast<char>(std::toupper(*(this->seq + pos + 2)));
            return (second_base != 'G') && (third_base == 'G');
        } else if (base == 'G') {
            if (pos <= 1) return 0;
            char second_base = static_cast<char>(std::toupper(*(this->seq + pos - 1)));
            char third_base = static_cast<char>(std::toupper(*(this->seq + pos - 2)));
            return (second_base != 'C') && (third_base == 'G');
        }
        return 0;
    }
    int is_chh(int ref_pos) const {
        int pos = ref_pos - this->start_pos;
        if (pos >= this->length) return 0;
        char base = static_cast<char>(std::toupper(*(this->seq + pos)));
        if (base == 'C') {
            if (pos + 2 >= this->length) return 0;
            char second_base = static_cast<char>(std::toupper(*(this->seq + pos + 1)));
            char third_base = static_cast<char>(std::toupper(*(this->seq + pos + 2)));
            return (second_base != 'G') && (third_base != 'G');
        } else if (base == 'G') {
            if (pos <= 1) return 0;
            char second_base = static_cast<char>(std::toupper(*(this->seq + pos - 1)));
            char third_base = static_cast<char>(std::toupper(*(this->seq + pos - 2)));
            return (second_base != 'C') && (third_base != 'C');
        }
        return 0;
    }
    int get_start_pos() const { return this->start_pos; }
    int get_length() const { return this->length; }
};

#endif // MODALIB_REFERENCE_H
