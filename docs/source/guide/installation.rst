Installation
============

Compiling HTSlib on Windows
---------------------------

.. code-block:: bash

    pacman -S mingw-w64-i686-libtool
    pacman -S mingw-w64-i686-toolchain
    pacman -S automake
    pacman -S autoconf

.. code-block:: bash

    autoreconf -i
    ./configure --disable-bz2 --disable-lzma --disable-dependency-tracking
    make
    make check
    make install


BAMConvertor
------------

You will need to put HTSlib in the ``dependencies`` folder. HTSlib can either be:

#. found at ``/lustre1/project/stg_00019/research/Antoine/MODaLiB/step1-bam-processing/dependencies/htslib-1.13/`` on the VSC, or
#. downloaded and unzipped at ``./dependencies/htslib-1.13/`` from `https://github.com/samtools/htslib/releases/ <https://github.com/samtools/htslib/releases/>`_.

First compile HTSlib:

.. code-block:: bash

    cd /path/to/dependencies/htslib-1.13/
    make clean
    make

Compile the project found in current folder:

.. code-block:: bash

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/dependencies/htslib-1.13/
    cd /lustre1/project/stg_00019/research/Antoine/MODaLiB/step1-bam-processing
    make clean
    make


Data fusion tool
----------------

First install Python dependencies:

.. code-block:: bash

    pip3 install -r requirements.txt


Compile the Cython and C++ code:

.. code-block:: bash

    cd step2-fragment-level-analysis
    python3 setup.py install --user

