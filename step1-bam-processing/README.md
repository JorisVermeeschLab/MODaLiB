Process BAM files and save fragment-level information
=====================================================

## Installation

You will need to put HTSlib in the `dependencies` folder. HTSlib can either be:
- found at `/lustre1/project/stg_00019/research/Antoine/MODaLiB/step1-bam-processing/dependencies/htslib-1.13/` on the VSC, or
- downloaded and unzipped at `./dependencies/htslib-1.13/` from [https://github.com/samtools/htslib/releases/](https://github.com/samtools/htslib/releases/).

First compile HTSlib:
```bash
cd /path/to/dependencies/htslib-1.13/
make clean
make
```

Compile the project found in current folder:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/dependencies/htslib-1.13/
cd /lustre1/project/stg_00019/research/Antoine/MODaLiB/step1-bam-processing
make clean
make
```

---

## Usage

```bash
BAMConvertor perRead path/to/reference/genome path/to/bam/file -o output/file --OT ot --OB ob'
```

The reference genome can be found at `/lustre1/project/stg_00019/methylation_data/human_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa`.
The inclusion bounds `ot` and `ob` can be used to mitigate the methylation biases at the ends of the reads.

The output file is written in a custom binary format described in the next section.

---

## Output file format

The file can be decomposed in chunks of 22 bytes. Each chunk follows the following format:
- Start chromosome (uint8_t, 1 byte)
- Start position (uint32_t, 4 bytes)
- End chromosome (uint8_t, 1 byte)
- End position (uint32_t, 4 bytes)
- Insert size (uint16_t, 2 bytes)
- Number of methylated CpG sites spanned by the fragment (uint16_t, 2 bytes)
- Number of unmethylated CpG sites spanned by the fragment (uint16_t, 2 bytes)
- Number of point differences with the reference genome (uint16_t, 2 bytes)
- Number of cytosines + guanines (uint16_t, 2 bytes)
- Number of informative bases (uint16_t, 2 bytes)
