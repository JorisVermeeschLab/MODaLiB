MODaLiB
=======

Multimodal analysis framework for cancer detection in cell-free DNA (cfDNA) samples.

Installation
------------

```bash
conda env create -f env.yml
conda activate modalib
pip install ensembl-converter qpsolvers[cvxopt] MLstatkit torch scanpy gseapy statsmodels
pip install snakemake pulp=2.7.0
conda install bioconda::methyldackel
```

Install GCparagon dependencies and clone the repository to any location:
[https://github.com/BGSpiegl/GCparagon](https://github.com/BGSpiegl/GCparagon)

Install HTSlib and add it to the default library path:
```bash
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:dependencies/htslib-1.18/
```

Compile the first part of the pipeline (`BAMConvertor`) with:
```bash
cd src/step1-bam-processing
make clean
make
cd ../..
```

Compile the Cython extensions with:
```
cd src/step2-fragment-level-analysis
python setup.py install --user
```

Build reference gene expression matrix in NumPy format
------------------------------------------------------

First, download The Tabula Sapiens dataset (all cell types) into the `gene-expression/` folder.
The dataset is available at:
[https://datasets.cellxgene.cziscience.com/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad](https://datasets.cellxgene.cziscience.com/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad).

Make sure that the `ensembl-ids.txt` file is in that directory.
It contains the mapping between Ensembl IDs and gene symbols, and can be obtained from:
[https://www.biotools.fr/human/ensembl_symbol_converter](https://www.biotools.fr/human/ensembl_symbol_converter).


Usage
-----

```bash
snakemake --cores 8 --config benchmark=true
```
