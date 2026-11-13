MODaLiB
=======

Multimodal analysis framework for cancer detection in cell-free DNA (cfDNA) samples.

Installation
------------

```bash
conda env create -f env.yml
conda activate modalib
pip install ensembl-converter qpsolvers[cvxopt] MLstatkit torch scanpy gseapy statsmodels Cython
pip install snakemake pulp==2.7.0
pip install plotly==5.24.1 kaleido==0.2.1
pip install finaletoolkit numpy==1.23
conda install bioconda::methyldackel
```

```bash
finaletoolkit gap-bed hg38 regions/hg38.gaps.bed
```

Install GCparagon dependencies and clone the repository to any location:
[https://github.com/BGSpiegl/GCparagon](https://github.com/BGSpiegl/GCparagon)


Compile the first part of the pipeline (`BAMConvertor`) with:
```bash
export HTSLIB=/lustre1/project/stg_00019/research/Antoine/MODaLiB/step1-bam-processing/dependencies/htslib-1.18/
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

Run the standard workflow with:

```bash
snakemake --cores 8
```

After configuring all contender commands as described below, run:

```bash
snakemake --cores 8 --config benchmark=true
```

Contender benchmarks
--------------------

Benchmark rules are not added to the workflow unless `benchmark=true`. When
enabled, FinaleToolkit runs DELFI and WPS directly. Set `benchmark_inputs` in
`config/config.yaml` to the chromosome sizes file, DELFI bin BED, and WPS site
BED. Install FinaleToolkit with `pip install finaletoolkit`.

Each argument can use `{sample}`, `{bam}`, `{reference}`, `{genome_2bits}`, and
`{output_dir}`. Configure `benchmark_commands` with argument lists for
`griffin`, `cancerdetector`, and `finaleme`; include each tool's model files and
other required options in its command. Commands run without a shell, once per
sample, on the GCparagon-corrected BAM. Each method gets its own output
directory under `results/benchmarks/<sample>/`. The workflow records Snakemake
resource/timing statistics and writes a JSON manifest only after all methods
finish successfully. Configure and install each tool before enabling the
benchmark flag.
