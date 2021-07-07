# Spatial transcriptomics sequencing

## Structure of the pipeline

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

## Inputs

### Samples

Each sample can be either defined in the `projects` variable in the `config.yaml` or in the `additional_projects`. The first expects a list of elements with each having two variables: `flowcell_id` and `sample_sheet`. The second expects several variables: `project_id`, `sample_id`, `puck_id`, `species`, `R1` and `R2`



## Configuration

The pipeline is implemented in snakemake. All metadata of the experiments (experiment\_name, flowcell\_id, species, etc) should be put in a `config.yaml` file. An example `config.yaml` file is in the root of this repo.

To run the snakemake script, the `snakemake` python library is required (installed with `pip` or `conda`). The script requires at least 6 threads to run, this is due to pipeing several commands one after the other to descrease runtime.

**Example run:**

`snakemake --snakefile path_to_snakefile --configfile path_to_configfile`.

This will create the output in the directory in which the command is run. Note, that all samplesheet-flowcell_id paris should be ideally in one configfile somewhere.

### Produced directory structure

The following directory structure will be produced by the snakemake file

        .
        |-- demultiplex_data                            # demultiplexed data folders, one per each samplesheet
        |   |-- 200110_STS_017_4-7-8STS_018_1           # directory names are identical to the samplesheet names
        |   |   |-- Stats
        |   |   |-- sts_017
        |   |   |-- sts_018
        |   |   |-- Undetermined_S0_R1_001.fastq.gz
        |   |   `-- Undetermined_S0_R2_001.fastq.gz
        |   `-- 20191206_spatseq_smples3-4
        |       |-- indicator.log
        |       |-- Reports
        |       |-- Stats
        |       |-- sts_0xxx
        |       |-- Undetermined_S0_R1_001.fastq.gz
        |       `-- Undetermined_S0_R2_001.fastq.gz
        |-- sts_017                                     # root output directory, one per project
        |   |-- data
        |   |   |-- sts_017_4                           # directory containing results of running the pipeline. one per sample 
        |   |   |-- sts_017_7
        |   |   `-- sts_017_8
        |   `-- reads                                   # reads directory, one per sample
        |       |-- fastqc
        |       |-- raw
        |       `-- reversed
        |-- sts_018
        |   |-- data
        |   |   `-- sts_018_1
        |   `-- reads
        |       |-- fastqc
        |       |-- raw
        |       `-- reversed
        `-- sts_0xxx
            |-- data
            |   |-- sts_01
            |   `-- sts_02
            `-- reads
                |-- fastqc
                |-- raw
                `-- reversed
