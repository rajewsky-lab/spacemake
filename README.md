# Spatial transcriptomics sequencing

This repository collects all scripts and tools used for analyzing the sequencing side of the spatial transcriptomics datasets. The following steps are currently performed:

## Demultiplex the data
This assumes that the sample sheet has been provided and that the raw data has been copied to the basecalls folder. The tool `bcl2fastq` is used to demultiplex the data.

## Rename the fastq files
It is important to rename the `.fastq` files so that the namings are meaningful.

## Reverse the fastq files
Read 1 needs to be reversed to match the barcodes of the optical side.

## Run FastQC on the fastq files
Run it on all files. Do QC.

## Run the sequencing analysis pipeline
After this the sequences are analyzed. It needs to be provided the (i) species to map onto and (ii) the filename of the sample. 

## Produce the QC sheet of the sequencing data
After everything is finished, a `python` script (containing an `R` part) is being run to produce the QC sheet for the sample. There's the `qc_sequencing_parameters.yaml` file which contains metadata for the experiment/sample and currently needs to be created automatically. Could be automized, with taking info partially from the sample sheet.
