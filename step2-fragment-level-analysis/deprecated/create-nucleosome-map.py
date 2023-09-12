# -*- coding: utf-8 -*-
#
#  create-nucleosome-map.py
#  
#  Copyright 2022 Antoine Passemiers <antoine.passemiers@gmail.com>
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

import os

from modalib.nucleosomes import save_piled_up_fragment_counts, create_reference_nucleosome_map

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(ROOT, 'step1-bam-processing/out')
folder = os.path.join(ROOT, '../../nucleosomes')


FILENAMES = [
    'GC119350',  'GC112864', 'GC119351', 'GC119352', 'GC112896', 'GC119353',
    'GC119354', 'GC119355', 'GC119356', 'GC110675', 'GC1003986', 'GC1003988',
    'GC112897', 'GC112898'
]
# FILENAMES = ['GC110680']


if not os.path.isdir(folder):
    os.makedirs(folder)

if not os.path.exists(os.path.join(folder, 'refcounts.npz')):
    input_filepaths = []
    for filename in FILENAMES:
        input_filepaths.append(os.path.join(DATA_PATH, filename))
    save_piled_up_fragment_counts(input_filepaths, os.path.join(folder, 'refcounts.npz'))
else:
    create_reference_nucleosome_map(os.path.join(folder, 'refcounts.npz'), folder)
