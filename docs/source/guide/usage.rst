Usage
=====

Process BAM files
-----------------

.. code-block:: bash

    BAMConvertor perRead path/to/reference/genome.fa path/to/file.bam -o output/file.modalib --OT ot --OB ob

The reference genome can be found at ``/lustre1/project/stg_00019/methylation_data/human_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa``.
The inclusion bounds ``ot`` and ``ob`` can be used to mitigate the methylation biases at the ends of the reads.

By default, reads are assumed to be paired. If not, the ``--single_read`` flag should be added.

By default, unmethylated cytosines are assumed to be converted to uracils (i.e., using bisulfite sequencing).
Alternatively, methylation information can be stored in the SAM/BAM optional fields called ``MM`` and ``ML``.
The specification of SAM/BAM auxiliary fields is available at 
`https://samtools.github.io/hts-specs/SAMtags.pdf <https://samtools.github.io/hts-specs/SAMtags.pdf>`_.
In order to use the auxiliary fields to call methylation, the ``--useMMML`` flag should be added to the command line.
Only the ``C+h`` base modification is considered: all other modification types will be ignored.
A CpG site will be considered as methylated if and only if the corresponding confidence value stored in the ``ML`` field
exceeds a given threshold. This threshold can be set with the optional command line argument ``--mlCutoff cutoff``.
``cutoff`` should be an integer between 1 and 255. Let's note that methylated cytosines at loci that do not correspond to a CpG site will be ignored.

.. CAUTION::

    Inclusion bounds are not yet compatible with the use of the ``--useMMML`` flag.


Modalib File format
-------------------

The output file is written in a custom binary format that is faster to process by the data fusion tool.
The file can be decomposed in chunks of 22 bytes. Each chunk follows the following format:

#. Start chromosome (uint8_t, 1 byte)
#. Start position (uint32_t, 4 bytes)
#. End chromosome (uint8_t, 1 byte)
#. End position (uint32_t, 4 bytes)
#. Insert size (uint16_t, 2 bytes)
#. Number of methylated CpG sites spanned by the fragment (uint16_t, 2 bytes)
#. Number of unmethylated CpG sites spanned by the fragment (uint16_t, 2 bytes)
#. Number of point differences with the reference genome (uint16_t, 2 bytes)
#. Number of cytosines + guanines (uint16_t, 2 bytes)
#. Number of informative bases (uint16_t, 2 bytes)

Convention: filenames extensions should be `.modalib`.


Reference nucleosome map
------------------------

TODO


Extract features
----------------

To perform multimodal analysis of cell-free DNA in the old fashion, genomic and epigenomic features should
be extracted from each sample separately:

.. code-block:: python

    from modalib.core import Sample, extract_features

    samples = [
        Sample('path/to/sample/1'),
        Sample('path/to/sample/2'),
        Sample('path/to/sample/3')
    ]
    matrices = extract_features(samples, ref_nuc_filepath)

`ref_nuc_filepath` is the path to the reference nucleosome map.
For example, `/lustre1/project/stg_00019/research/Antoine/MODaLiB/step2-fragment-level-analysis/nucleosomes/Homo_sapiens.hsNuc0100101.nucleosome.DANPOSPeak.sorted.bed`.

The output ``matrices`` is a `dict` with as many entries as chromosomes. For example, the key for chromosome 6 is simply 6.
Each value is a list of size `p`, where `p` is the number of feature types (e.g., read counts, methylation).
Each element of such list is itself a list of size `n`, where `n` is the number of input samples ``len(samples)``.
Finally, each element of such sub-list is an array of size `m`, where `m` is the number of bins in corresponding chromosome (i.e., of size 50 kb).

The 5 feature types are, by order of appearance:

#. Average nucleosome positioning score
#. Ratio of the number of long fragments (>150 bp) over short fragments
#. Methylation ratio
#. Fragment count
#. "Mutation" burden (fraction of bases differing from the reference genome)


Compute fragment-level kernel matrices
--------------------------------------

To perform multimodal analysis of cfDNA using fragment-level kernel matrices instead, the python code is very similar:

.. code-block:: python

    from modalib.core import Sample, pairwise_kernels

    samples = [
        Sample('path/to/sample/1'),
        Sample('path/to/sample/2'),
        Sample('path/to/sample/3')
    ]
    matrices = pairwise_kernels(filepaths, ref_nuc_filepath)

The output ``matrices`` is a `dict` with as many entries as chromosomes. For example, the key for chromosome 6 is simply 6.
Each value is a list of size `p`, where `p` is the number of kernels (e.g., read counts, methylation), which is the same as the number
of feature types listed in previous section.
Each element of such list is a matrix of size `n × n`, where `n` is the number of input samples ``len(samples)``.

The kernels are listed in the same way as the feature types listed in previous section. For example, the third kernel matrix from each list
encodes the sample-to-sample similarity of methylation patterns.
