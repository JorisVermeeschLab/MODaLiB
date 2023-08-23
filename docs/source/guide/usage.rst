Usage
=====

Process BAM files
-----------------

.. code-block:: bash

    BAMConvertor perRead path/to/reference/genome.fa path/to/file.bam -o output/file.modalib --OT ot --OB ob

The reference genome can be found at ``/lustre1/project/stg_00019/methylation_data/human_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa``.
The inclusion bounds ``ot`` and ``ob`` can be used to mitigate the methylation biases at the ends of the reads.

By default, reads are assumed to be paired. If not, the ``--singleRead`` flag should be added.

By default, unmethylated cytosines are assumed to be converted to uracils (i.e., using bisulfite sequencing).
Alternatively, methylation information can be stored in the SAM/BAM optional fields called ``MM`` and ``ML``.
The specification of SAM/BAM auxiliary fields is available at 
`https://samtools.github.io/hts-specs/SAMtags.pdf <https://samtools.github.io/hts-specs/SAMtags.pdf>`_.
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

Convention: filenames extensions should be ``.modalib``.


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


.. code-block:: python

    from modalib.nucleosomes import create_reference_nucleosome_map

    create_reference_nucleosome_map('folder/refcounts.npz', 'output/folder')


Extract features
----------------

.. code-block:: python

    from modalib.core import Sample, extract_features

    samples = [
        Sample('path/to/sample/1'),
        Sample('path/to/sample/2'),
        Sample('path/to/sample/3')
    ]
    nuc_folder = '/lustre1/project/stg_00019/research/Antoine/MODaLiB/step2-fragment-level-analysis/nucleosomes'
    matrices = extract_features(samples, nuc_folder)

``nuc_folder`` is a folder 

Compute fragment-level kernel matrices
--------------------------------------

.. code-block:: python

    from modalib.core import Sample, pairwise_kernels

    TODO
