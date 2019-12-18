# Spatial transcriptomics sequencing

This repository collects all scripts and tools used for analyzing the sequencing side of the spatial transcriptomics datasets. The following steps are currently performed:

## Demultiplex the data
This assumes that the sample sheet has been provided and that the raw data has been copied to the basecalls folder. The tool `bcl2fastq` is used to demultiplex the data.

## Rename the fastq files
It is important to rename the `.fastq` files so that the namings are meaningful. This is currently done manually.

## Reverse the fastq files
Read 1 needs to be reversed to match the barcodes of the optical side. Sometimes we sequence 20 bases and sometimes 21. Currently the script is hard-coded to select for the appropriate length, needs to be automized. We should also at this step complement the barcodes to match with the optical side without doing further manipulations there.

## Run FastQC on the fastq files
Run it on all files and spit the files. Do QC.

## Run the dropseq toolkit
After this the standard dropseq pipeline is being run. It needs to be provided the (i) species to map onto, (ii) the filename of the sample and (iii) the number of beads expected. The last is sometimes creating errors when there's not a clear knee plot in the data. Needs to be better automized.

## Produce the QC sheet of the sequencing data
After everything is finished, a `python` script (containing an `R` part) is being run to produce the QC sheet for the sample. There's the `parameters.yaml` file which contains metadata for the experiment/sample and currently needs to be created automatically. Could be automized.
