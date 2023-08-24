# -*- coding: utf-8 -*-
# setup.py

import os
import sys
import glob
from distutils.core import setup
from distutils.extension import Extension
try:
    from Cython.Build import cythonize
except ImportError:
    import Cython.Build
    cythonize = Cython.Build.cythonize
import numpy as np
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup as np_setup


# Project-dependant variable
PACKAGE_NAME = 'modalib'
SOURCE_FOLDER = 'modalib'
LANGUAGE_LEVEL = '3'  # Python 3
INCLUDE_DIRS = ['include']

from sys import platform
if platform == 'win32':
    EXTRA_COMPILE_ARGS = [
        '/std:c++20',
        '/Wall',
        '/Ox',
        '/fp:fast',
        '/openmp'
    ]
else:
    EXTRA_COMPILE_ARGS = [
        '-std=c++20',
        '-Wall',
        '-Wextra',
        '-O3',
        # '-O0', "-ggdb3",  # gdb --ex r --args python3 file.py <args>
        '-ffast-math',
        '-march=native',
        # '-fopenmp'
    ]

# Allowed source file extensions
ALLOWED_EXTENSIONS = {'.py', '.c', '.cpp'}

# Look for Cython source files
source_filepaths = [match for filepath in os.walk(SOURCE_FOLDER) \
        for match in glob.glob(os.path.join(filepath[0], '*.pyx'))]

# Create Cython extensions (only for generating the C code)
blacklist = list()
extensions = list()
for source_filepath in source_filepaths:
    sources = [source_filepath]
    extension_name = os.path.splitext(source_filepath)[0]
    extension_name = extension_name.replace('/', '.').replace('\\', '.')
    print('[CYTHON EXTENSION] %s %s' % (str(extension_name), str(sources)))
    extensions.append(
        Extension(
            extension_name,
            sources,
            include_dirs=[np.get_include()] + INCLUDE_DIRS))
    blacklist.append(os.path.splitext(source_filepath)[0] + '.c')
    blacklist.append(os.path.splitext(source_filepath)[0] + '.cpp')

GOT_BUILD_CMD = (sys.argv[1] in {'install', 'build', 'build_ext'})

if GOT_BUILD_CMD:

    # Remove __init__.py file if present in root folder.
    init_path = os.path.join(os.path.realpath(__file__), '../__init__.py')
    if os.path.isfile(init_path):
        os.remove(init_path)
        print('[NOTE] Root __init__.py file removed')

    # Generates the C files, but does not compile them
    # https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#cythonize-arguments
    compiler_directives = {'embedsignature': True, 'language_level' : LANGUAGE_LEVEL}
    cythonize(extensions, compiler_directives=compiler_directives, force=True)

# Look for sub-packages
start = os.path.join(SOURCE_FOLDER, '.')
sub_packages = set()
for sub_package in set([os.path.relpath(x[0], start=start).replace('\\', '.').replace('/', '.') for x in os.walk(SOURCE_FOLDER)]):
    if '__pycache__' not in sub_package and sub_package != '.':
        print('[SUB-PACKAGE] %s' % str(sub_package))
        sub_packages.add(PACKAGE_NAME + '.' + sub_package)

# Look for source files inside sub-packages
non_cython_source_files = list()
whitelist = list()
source_files = list()
for x in os.walk(SOURCE_FOLDER):
    sub_folder = os.path.relpath(x[0], start=start)
    for filename in x[2]:
        if sub_folder != '.':
            rel_filepath = os.path.join(sub_folder, filename)
        else:
            rel_filepath = filename
        if os.path.splitext(rel_filepath)[1] in ALLOWED_EXTENSIONS and '__init__' not in rel_filepath:
            module_name = os.path.splitext(rel_filepath)[0].replace('\\', '.').replace('/', '.')
            if not 'setup' in module_name:
                whitelist.append(os.path.join(SOURCE_FOLDER, rel_filepath))
                if whitelist[-1] in blacklist:
                    source_files.append(([rel_filepath], module_name))
                    print('[EXTENSION] %s' % str(source_files[-1]))

# List all the source files that are not a byproduct of cythonization
non_cython_source_files = list(set(whitelist) - set(blacklist))

# Compile and create package
config = Configuration(
        package_name=PACKAGE_NAME,
        parent_name=None,
        package_path=SOURCE_FOLDER,
        top_path=SOURCE_FOLDER)
for sub_package in sub_packages:
    config.add_subpackage(sub_package)

# Create one extension per cython file (not per C/C++ file),
# since python cannot call C/C++ code.
# For each extension, all sources that are not the byproduct
# of cythonization are added.
libraries = ['m'] if os.name == 'posix' else list()
include_dirs = INCLUDE_DIRS + [np.get_include()]
for sources, extension_name in source_files:
    sources = [os.path.join(SOURCE_FOLDER, source) for source in sources]
    sources += non_cython_source_files
    sources = list(set(sources))
    print('SOURCES', extension_name, sources)
    config.add_extension(
        extension_name,
        sources=sources,
        include_dirs=include_dirs+[SOURCE_FOLDER],
        libraries=libraries,
        extra_compile_args=EXTRA_COMPILE_ARGS)

if GOT_BUILD_CMD:
    # os.environ['CFLAGS'] = ' '.join(EXTRA_COMPILE_ARGS)
    np_setup(**config.todict())
