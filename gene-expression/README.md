Build reference gene expression matrix in NumPy format
======================================================

First, download The Tabula Sapiens dataset (all cell types) at:
[https://datasets.cellxgene.cziscience.com/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad](https://datasets.cellxgene.cziscience.com/10df7690-6d10-4029-a47e-0f071bb2df83.h5ad).

Make sure that the `ensembl-ids.txt` file is in the working directory.
It contains the mapping between Ensembl IDs and gene symbols, and can be obtained from:
[https://www.biotools.fr/human/ensembl_symbol_converter](https://www.biotools.fr/human/ensembl_symbol_converter).

Next, run all the code cells in the `tabula-sapiens-to-numpy.ipynb` notebook to produce the `expression-profiles.npz` file.
