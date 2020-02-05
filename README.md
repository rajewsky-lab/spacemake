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
After everything is finished, a `python` script (containing an `R` part) is being run to produce the QC sheet for the sample. There's the `qc_sequencing_parameters.yaml` file which contains metadata for the experiment/sample and currently needs to be created automatically. Could be automized, with taking info partially from the sample sheet.

## Snakemake

The pipeline is implemented in snakemake. All metadata of the experiments (experiment\_name, flowcell\_id, species, etc) should be put in a `config.yaml` file. An example `config.yaml` file is in the root of this repo.

To run the snakemake script, the `snakemake` python library is required (installed with `pip` or `conda`). The script requires at least 6 threads to run, this is due to pipeing several commands one after the other to descrease runtime.

**Example run:**

`snakemake --snakefile path_to_snakefile --configfile path_to_configfile`.

This will create the output in the directory in which the command is run. Note, that all samplesheet-flowcell_id paris should be ideally in one configfile somewhere.

### Produced directory structure

The following directory structure will be produced by the snakemake file

    .
    └── <run_name_1>
        ├── data
        │   ├── <experiment_name_1>
        │   │   ├── dge             # folder containing all DGEs
        │   │   ├── qc_sheet        # folder with the qc sheet
        │   │   └── reports         # folder with all report, and summary files from the pipeline
        │   └── <experiment_name_2>
        │       ├── dge
        │       ├── qc_sheet
        │       └── reports
        ├── demux_data              # demultiplexing root directory
        │   ├── Reports
        │   │   └── html
        │   ├── Stats
        │   └── <project_name_1>
        │       ├── sts_01
        │       └── sts_02
        └── reads                   # reads root directory
            ├── fastqc              # fastqc folder
            ├── raw                 # directory containing symbolic links to the demultiplexed reads
            └── reversed            # directory containing reversed R1 and symbolic link to R2 from raw reads

