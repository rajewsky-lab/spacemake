import importlib.metadata

__version__ = "0.9.1"
__author__ = [
    "Nikos Karaiskos",
    "Tamas Ryszard Sztanka-Toth",
    "Marvin Jens",
    "Daniel Leon-Perinan",
]
__license__ = "GPL"
__email__ = [
    "nikolaos.karaiskos@mdc-berlin.de",
    "tamasryszard.sztanka-toth@mdc-berlin.de",
    "marvin.jens@charite.de",
    "daniel.leonperinan@mdc-berlin.de",
]

author_contributions = """
Spacemake is built on snakemake scripts originally developed by Nikos Karaiskos
for the analysis of dropseq data. These gradually evolved into a robust workflow for
spatial transcriptomics data analysis that was improved and generalized to work
with different ST technologies by Tamas Ryszard Sztanka-Toth. Marvin Jens contributed
longread analysis code and support for converting fastq to BAM as a first step.
Many features of the automated analysis and integration with Novosparc were added by
Tamas, in close collaboration with Nikos, culminating in the first spacemake
publication:

    https://doi.org/10.1093/gigascience/giac064

Marvin then added new building blocks to successively replace the java-based 
dropseq tools with python/pysam based code: cutadapt_bam.py, annotator.py, as well
as the ability to align raw reads to multiple indices, in close collaboration
with Nikos & Tamas.

Spacemake is actively maintained by Dani, Marvin and Nikos.
"""

roadmap = [
    (
        "0.5.5",
        "universal ST support and utility, "
        "novosparc integration. Sztanka-Toth et al. 2022",
    ),
    (
        "0.7",
        "support multiple mapping indices, "
        "bulk samples, custom user-defined snakemake rules",
    ),
    ("0.7.7", "much lower RAM usage and speed boost"),
    (
        "0.8",
        "much faster cmdline. 10x speed-up for fastq_to_uBAM. "
        "Lower memory for BamTagHistogram replacement",
    ),
    (
        "0.9",
        "CRAM on-disk format together with other tweaks "
        "-> disk footprint down to 1/3 of 0.8. "
        "Requires migration of directory structure via "
        "new `spacemake migrate` command."
        "Replace DropSeq pre-processing tools with cutadapt "
        "code for polyA-trimming and adapter removal.",
    ),
    (
        "0.9.1",
        "QC reports are now produced with python/notebooks. "
        "Completely dropped dependencies on all R packages."
    ),
    (
        "0.9.2",
        "also count alignments against non-genome indices, "
        "using scbamtools.count. Spatial output is merged "
        "from genome and non-genome alignments.",
    ),
    (
        "1.x",
        "completely replace dropseq tools. "
        "Own annotator and towards entirely scanpy workflow.",
    ),
    ("1.x", "efficient handling of 1E8+ spatial barcodes (seq-scope etc.)"),
    ("1.x", "add interactive data exploration support (shiny?)"),
    ("2.x", "cmdline interface cleanup and remote API support"),
    ("2.x", "cython magic to speed up parallel BAM processing via shared memory"),
]
