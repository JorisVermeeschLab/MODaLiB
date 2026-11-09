//  context.h
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

#ifndef MODALIB_CONTEXT_H
#define MODALIB_CONTEXT_H

#include <stdint.h>


enum Context {
    NOT_C_CONTEXT,
    CpG_CONTEXT,
    CHG_CONTEXT,
    CHH_CONTEXT,
    UNKNOWN_CONTEXT
};


bool is_cytosine(Context context);

uint8_t as_uint8_t(Context context);


#endif  // MODALIB_CONTEXT_H
