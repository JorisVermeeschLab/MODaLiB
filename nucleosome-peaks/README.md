How to get the reference nucleosome map
=======================================

First download the `GSE71378_CA01.bb` file available at
[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71378),
and place it in the current directory.

Next, extract the nucleosome core centers with:

```bash
python extract-peaks.py
```

Finally, convert the genomic coordinates to hg38 with:

```bash
python convert-peaks-to-hg38.py
```
