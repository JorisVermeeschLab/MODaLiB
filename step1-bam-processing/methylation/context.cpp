//  context.cpp
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

#include "context.h"


bool is_cytosine(Context context) {
    switch (context) {
        case NOT_C_CONTEXT:
            return false;
        case CpG_CONTEXT:
        case CHG_CONTEXT:
        case CHH_CONTEXT:
            return true;
        default:
            return false;
    }
}


uint8_t as_uint8_t(Context context) {
    switch (context) {
        case NOT_C_CONTEXT:
            return 0;
        case CpG_CONTEXT:
            return 1;
        case CHG_CONTEXT:
            return 2;
        case CHH_CONTEXT:
            return 3;
        default:
            return 0;
    }
}
