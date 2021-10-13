# spacemake: pipeline for processing and analysing sequencing based spatial-transcriptomics data

## Installation

To install spacemake type `pip install spacemake`

This repository collects all scripts and tools used for analyzing the sequencing side of the spatial transcriptomics datasets. The following steps are currently performed:

### Demultiplex the data
This assumes that the sample sheet has been provided and that the raw data has been copied to the basecalls folder. The tool `bcl2fastq` is used to demultiplex the data.

### Extract Bead Barcode and UMI from R1 and R2
Depending on the R1 and R2 structure, which can be defined in the top-level config.yaml, spacemake will create a tagged but unmapped.bam, with CB (Bead Barcode) and MI (UMI) tags. This .bam file will then be used as an input by STAR

### Processing

Once we extracted the CB and MI tags, the data is processed as follows:
- Adapter removal (R2)
- PolyA removal (R2)
- Map reads with STAR (species can be defined in either the samplesheet or in the config.yaml)
- Postprocessing
- Join data with the bead position information
- Analyse the data automatically, produce reports
