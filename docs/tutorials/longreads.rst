Longreads Integration
=====================

A Spacemake sample can not only be associated with Illumina sequencing reads 
via the ``--R1`` and ``--R2`` command line arguments, but can also be assigned 
long reads (e.g. PacBio, Nanopore, etc.) using the ``--longreads`` command-line 
argument.

We have added this functionality for the purpose of trouble-shooting 
problems with library construction (see the Spacemake paper). Long sequencing reads can
capture the entire cDNA sequence, including Illumina sequencing adapters, primer handles, 
cell and UMI barcodes, polyA/T as well as the actual cDNA insert.

If you add long reads to a sample, you will enable spacemake to annotate a catalog of 
oligo sequences, thought of as building blocks of your library, against every long read. 
This allows to assess the following feaures:

  - the fraction of molecules which conform to the expected layout of building blocks.
    For example, SMART-handle, (barcode), polyT, (cDNA insert), TSO-handle. Parts in 
    parenthesis here are inferred from the spacing between the adjacent blocks but not 
    directly annotated

  - the fraction of molecules which are missing one or more the expected building blocks
  
  - size distributions of all matches and distributions of their start and end positions,
    which in turn allows to infer the size distributions of inserts and barcode sequences.

  - *concatamerizations* and unexpected, *multiple* occurrences of any building block, 
  - pointing to undesired side-reactions.


Building blocks and Signatures
------------------------------

In order to fully utilize the longread functionality, spacemake needs to know about the sequences 
of all building blocks, as well as their expected layout. For the building blocks, 
spacemake searches a FASTA formatted file in 
``<spacemake-root>/spacemake/data/oligo_blocks.fa``. 
Here is an excerpt from this file:

.. code-block::

   >P5
   AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGG

   >SMART_Primer
   AAGCAGTGGTATCAACGCAGAGT

   >bead_start
   AAGCAGTGGTATCAACGCAGAGTAC

   >dN-SMRT
   AAGCAGTGGTATCAACGCAGAGTGA

   >TSO
   AAGCAGTGGTATCAACGCAGAGTGAATGGG

   >polyT
   TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT


You can see here that ``SMART_Primer`` and ``bead_start`` are very similar sequences, 
differing only at the end, with ``bead_start`` containing an extra ``AC`` dinucleotide. 
In case of multiple (partial) matches to a building block, a match with higher score 
(ie more matches) will supersede an overlapping match with lower score. This allows to
distinguish the handle at the start of the capture oligos (which features the ``AC`` at the 3' end) 
from other occurrences of the SMART primer site, which might be introduced via a different route.

Tying it together is a ``signature``. Signatures are defined (at the moment manually) in 
a ``longread.yaml`` file. If spacemake does not find a file with that name at the root of 
your ``projects`` folder, it will default to loading ``<spacemake-root>/spacemake/data/longread.yaml``.

Here is an excerpt from this file describing the expected layout of dropseq/Nadia cDNA libraries:

.. code-block:: yaml

   signatures:
       dropseq:
           prio: -1
           label: "dropseq"
           read1_primer: bead_start
           read2_primer: N70X
           CB: "r1[8:20]"
           UMI: "r1[0:8]"
           color: "gray"
           intact_bead: "P5,bead_start,polyT,N70X"
           bead_related: "bead_start"
           ignore_matches: ["chr_bead_start", "chr_TSO","chrV2_RT_PRIMER","chrV3_RT_PRIMER", "OP1", "OP2", "OP2_2s", "OP3"]

Other pre-defined signatures at the moment are ``chromium`` and (the almost identical) ``visium``.

The ``prio`` and ``color`` fields are only relevant for overview plots across multiple samples 
and will affect the ordering (prio) and color for the visual representation of sample metrics.
Most important are the ``intact_bead`` and ``bead_related`` fields. 

The ``intact_bead`` field lists all buidling blocks expected to be present (in that order) 
on a correct library molecule. P5 and N70X are always considered optional at this point, 
because you may choose to perform long read sequencing on a library before (absent) or 
after index PCR (present).

The ``bead_related`` field specifies a single building block that you expect to be present in all
molecules that can reasonably be expected to derive from the used capture technology.
In the case of dropseq/nadia beads, this would be the SMART handle, followed by ``AC``. 
Occurrences of the ``bead_related`` building block are further used to orient every long read, 
as reads can be either in forward, or reverse-complement orientation.

What happens?
-------------

As soon as at least one of your samples is associated with long sequencing reads, ``spacemake run`` 
will invoke some dedicated tools. Specifically

   1. long reads will be *aligned* against all known building blocks
   2. (overlapping) matches will be *collected* and integrated for each read
   3. based on the presence/absence of each block, each read will be *classified*
   4. *statistics* on the observed long read classes will be gathered, with particular emphasis on the 
      reads falling into the class defined as ``intact_bead`` .
   5. *cDNA* will be extracted and mapped to the genome via `STAR-long`
   6. association of mappability and building block presence/absence is investigated
   7. report *plots* are generated for each sample in ``/processed_data/{sample_id}/pacbio/reports``

After these steps are completed for every sample with long reads, *overview plots* are generated, which 
present high level results across all samples, side-by-side in ``<projects folder>/pacbio_overview``.

If you like to utilize any of these functions outside of the spacemake/snakemake workflow you can either 
invoke the longread command via ``python -m spacemake.longread`` or by importing the ``spacemake.longread``
module from your own python scripts.

Example
-------

Here is a full example using a small test data-set. We will download the test data, add the 
sample to a spacemake project, run the analysis, and have a look at the output generated. We use ``wget`` 
to download a small test data set. Alternatively, you can use 
``git clone https://github.com/rajewsky-lab/spacemake-test-data.git`` to check out a collection of different
test data. Here, we further assume that a species labelled ``human`` has already been set up for spacemake [TODO: add link].

.. code-block::

   wget https://github.com/rajewsky-lab/spacemake-test-data/raw/main/longread/SRR9008425_subsample.fastq.gz

   spacemake projects add_sample --project_id=test --sample_id=test_longread --longreads=SRR9008425_subsample.fastq.gz --longread-signature=chromium --species=human
   spacemake run -np --cores=8

Spacemake lays out the following tasks to process our longread sample (shortened for brevity):

.. code-block::

    Job counts:
            count   jobs
            1       all
            1       cmd_align
            1       cmd_alnstats
            1       cmd_annotate
            1       cmd_edits
            1       cmd_extract
            1       cmd_overview
            1       cmd_report
            1       map_cDNA
            9
    This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.

    rule cmd_align:
        python -m spacemake.longread --parallel=8 --config=longread.yaml --cache=projects/test/processed_data/test_longread/pacbio/cache/ --annotation-out=projects/test/processed_data/test_longread/pacbio/annotation/ --stats-out=projects/test/processed_data/test_longread/pacbio/stats/ --report-out=projects/test/processed_data/test_longread/pacbio/reports/ --examples-out=projects/test/processed_data/test_longread/pacbio/examples/ --sample=test_longread --signature=chromium align SRR9008425_subsample.fastq.gz 

    rule cmd_extract:
        python -m spacemake.longread --parallel=1 --config=longread.yaml --cache=projects/test/processed_data/test_longread/pacbio/cache/     --annotation-out=projects/test/processed_data/test_longread/pacbio/annotation/     --stats-out=projects/test/processed_data/test_longread/pacbio/stats/     --report-out=projects/test/processed_data/test_longread/pacbio/reports/     --examples-out=projects/test/processed_data/test_longread/pacbio/examples/ --sample=test_longread --signature=chromium extract SRR9008425_subsample.fastq.gz 2> projects/test/processed_data/test_longread/pacbio/cDNA/test_longread.log > projects/test/processed_data/test_longread/pacbio/cDNA/test_longread.fa

    rule cmd_annotate:
        python -m spacemake.longread --parallel=1 --config=longread.yaml --cache=projects/test/processed_data/test_longread/pacbio/cache/     --annotation-out=projects/test/processed_data/test_longread/pacbio/annotation/     --stats-out=projects/test/processed_data/test_longread/pacbio/stats/     --report-out=projects/test/processed_data/test_longread/pacbio/reports/     --examples-out=projects/test/processed_data/test_longread/pacbio/examples/ --sample=test_longread --signature=chromium annotate SRR9008425_subsample.fastq.gz
    
    rule map_cDNA:
        mkdir -p projects/test/processed_data/test_longread/pacbio/cDNA/tmp/
        STARlong             --runThreadN 8             --genomeDir  species_data/human/star_index             --genomeLoad NoSharedMemory             --readFilesIn projects/test/processed_data/test_longread/pacbio/cDNA/test_longread.fa             --readFilesType Fastx             --outSAMtype BAM Unsorted             --outSAMunmapped Within             --outSAMattributes All             --outSAMprimaryFlag AllBestScore             --outStd BAM_Unsorted --outFilterMultimapScoreRange 2 --outFilterScoreMin 0 --outFilterScoreMinOverLread 0             --outFilterMatchNminOverLread 0             --outFilterMatchNmin 30             --outFilterMismatchNmax 1000             --winAnchorMultimapNmax 200             --seedSearchStartLmax 12             --seedPerReadNmax 100000             --seedPerWindowNmax 100             --alignTranscriptsPerReadNmax 100000             --alignTranscriptsPerWindowNmax 10000             --outFileNamePrefix projects/test/processed_data/test_longread/pacbio/cDNA/tmp |             /data/rajewsky/shared_bins/Drop-seq_tools-2.4.0/TagReadWithGeneFunction             I=/dev/stdin             O=projects/test/processed_data/test_longread/pacbio/cDNA/test_longread.bam             ANNOTATIONS_FILE=species_data/human/annotation.gtf

    rule cmd_report:
        python -m spacemake.longread     --parallel=1     --config=longread.yaml          --cache=projects/test/processed_data/test_longread/pacbio/cache/     --annotation-out=projects/test/processed_data/test_longread/pacbio/annotation/     --stats-out=projects/test/processed_data/test_longread/pacbio/stats/     --report-out=projects/test/processed_data/test_longread/pacbio/reports/     --examples-out=projects/test/processed_data/test_longread/pacbio/examples/ --sample=test_longread --signature=chromium report

    rule cmd_edits:
        python -m spacemake.longread     --parallel=1     --config=longread.yaml          --cache=projects/test/processed_data/test_longread/pacbio/cache/     --annotation-out=projects/test/processed_data/test_longread/pacbio/annotation/     --stats-out=projects/test/processed_data/test_longread/pacbio/stats/     --report-out=projects/test/processed_data/test_longread/pacbio/reports/     --examples-out=projects/test/processed_data/test_longread/pacbio/examples/     --sample=test_longread     --signature=chromium       edits SRR9008425_subsample.fastq.gz

    rule cmd_overview:
        python -m spacemake.longread     --parallel=1     --config=longread.yaml       overview --output pacbio_overview/ projects/test/processed_data/test_longread/pacbio/stats/test_longread.report.tsv


Ok, sounds good. Let's do this.