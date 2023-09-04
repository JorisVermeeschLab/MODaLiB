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


Create a reference nucleosome map
---------------------------------

As a first step, the fragments stored in the different files should be aggregated by piling them up.
The output is a NumPy `.npz` file with multiple arrays, each array corresponding to one chromosome.
Each array has the same size as the chromosome (in bp), and each value corresponds to the number of fragments spanning that locus.

.. code-block:: python

    from modalib.nucleosomes import save_piled_up_fragments_counts

    input_filepaths = [
        'sample1.modalib',
        'sample2.modalib',
        ...
        'sample135.modalib'
    ]
    save_piled_up_fragment_counts(input_filepaths, 'folder/refcounts.npz')

It is recommended using controls or samples that are representative of the same population.

After computing the fragment counts, the nucleosome positions are called by executing the following function:

.. code-block:: python

    from modalib.nucleosomes import create_reference_nucleosome_map

    create_reference_nucleosome_map('folder/refcounts.npz', 'nucleosomes/folder')

This function will store one `.npy` file per chromosome. For example, `6.npy` contains a NumPy array
with the inferred positions for the nucleosomes in chromosome 6.
This array can be loaded as such:

.. code-block:: python

    positions = np.load('nucleosomes/folder/6.npy')


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
    matrices = extract_features(samples, nuc_folder)

``nuc_folder`` is the path to the folder containing the nucleosome reference map described in previous section.
On the VSC, there is a reference map already built for paired-end EMseq samples, located at
``/lustre1/project/stg_00019/research/Antoine/MODaLiB/step2-fragment-level-analysis/nucleosomes``.

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
    matrices = pairwise_kernels(filepaths, nuc_folder)

The output ``matrices`` is a `dict` with as many entries as chromosomes. For example, the key for chromosome 6 is simply 6.
Each value is a list of size `p`, where `p` is the number of kernels (e.g., read counts, methylation), which is the same as the number
of feature types listed in previous section.
Each element of such list is a matrix of size `n × n`, where `n` is the number of input samples ``len(samples)``.

The kernels are listed in the same way as the feature types listed in previous section. For example, the third kernel matrix from each list
encodes the sample-to-sample similarity of methylation patterns.
