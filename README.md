# Spatial transcriptomics sequencing

## Structure of the pipeline

This repository collects all scripts and tools used for analyzing the sequencing side of the spatial transcriptomics datasets. The following steps are currently performed:

### Demultiplex the data
This assumes that the sample sheet has been provided and that the raw data has been copied to the basecalls folder. The tool `bcl2fastq` is used to demultiplex the data.

### Rename the fastq files
It is important to rename the `.fastq` files so that the namings are meaningful.

### Reverse the fastq files
Read 1 needs to be reversed to match the barcodes of the optical side.

### Run FastQC on the fastq files
Run it on all files. Do QC.

### Run the sequencing analysis pipeline
After this the sequences are analyzed. It needs to be provided the (i) species to map onto and (ii) the filename of the sample. 

### Produce the QC sheet of the sequencing data
After everything is finished, a `python` script is being run to produce the QC sheet for the sample. There's the `qc_sequencing_parameters.yaml` file which contains metadata for the experiment/sample and currently needs to be created automatically. Could be automized, with taking info partially from the sample sheet.

## Snakemake

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
