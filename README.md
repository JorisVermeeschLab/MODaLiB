MODaLiB
=======

Multimodal analysis framework for cancer detection in cell-free DNA (cfDNA) samples.

Installation
------------

Install GCparagon dependencies and clone the repository to any location:
[https://github.com/BGSpiegl/GCparagon](https://github.com/BGSpiegl/GCparagon)

Install HTSlib and add it to the default library path:
```bash
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:dependencies/htslib-1.18/
```

Install MODaLiB's dependencies with:
```bash
pip install -r requirements.txt
```

Compile the first part of the pipeline (`BAMConvertor`) with:
```bash
cd step1-bam-processing
make clean
make
cd ..
```

Compile the Cython extensions of the second part of the pipeline with:
```
cd step2-fragment-level-analysis
python setup.py install --user
```

Regarding the `modalibpy` package, it does not require any compilation.

Usage
-----


