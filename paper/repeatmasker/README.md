Estimate the repetitiveness of DMRs
===================================

MethMarkerDB DMRs have already been annotated in:
- `regions/metdecode/dmrs-methmarkerdb-colorectal-with-rmo.csv`
- `regions/metdecode/dmrs-methmarkerdb-breast-with-rmo.csv`

Last column of each file is the repetitiveness, computed as the degree of overlap between a DMR and the RepeatMasker mask.

In case you wish to re-generate these files from scratch, please follow the instructions below.

First, download the `hg38_2020_rmsk.bed` file from the UCSC Repeat Browser website at:
[https://repeatbrowser.ucsc.edu/data/](https://repeatbrowser.ucsc.edu/data/) (RepeatMasker annotations (bed files for human genome assemblies) hg38).
Exact download link is: [https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38/hg38_2020_rmsk.bed](https://hgdownload.soe.ucsc.edu/hubs/RepeatBrowser2020/hg38/hg38_2020_rmsk.bed). Place the file in the current directory.

Next, run the `compute-overlaps.py` script.
