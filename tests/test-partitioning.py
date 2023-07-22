# -*- coding: utf-8 -*-
#
#  test-partitioning.py
#  
#  Copyright 2023 Antoine Passemiers <antoine.passemiers@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

from modalib.core import coords_to_partition_


def test_compute_partition_id():

    assert coords_to_partition_(1, 8392844, 3) == 0
    assert coords_to_partition_(1, 1938402488, 3) == 0
    assert coords_to_partition_(12, 893443, 3) == 1
    assert coords_to_partition_(22, 930274, 3) == 2

