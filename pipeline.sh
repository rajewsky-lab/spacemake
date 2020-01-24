#!/bin/bash
#  
# Master pipeline script for the computational analysis of Spatial Transcriptomics Sequencing data
#
# author : Nikos Karaiskos
# email : nikolaos.karaiskos@mdc-berlin.de
#
# version 0.1.6
#
###############################################################################



###############################################################################
#
# Variable holding the folder where the toolkit is

toolkit_folder='/home/nkarais/src/git/sts-sequencing'

#
###############################################################################




###############################################################################
#
# Variables passed as arguments
#
# The sample sheet for demultiplexing. Needs to be a full path.
sample_sheet_full_location=$1

# From this the folder of the demultiplexed data is inferred
sample_sheet_filename=$(echo $sample_sheet_full_location | sed 's:.*/::')
folder=${sample_sheet_full_location/$sample_sheet_filename/}demultiplexed_data

# We also read the flowcell_id that was used for sequencing
flowcell_id=$(cat ${sample_sheet_full_location/$sample_sheet_filename/}flowcell_id.txt)

# The genome species that the reads will be mapped onto
species=$2

#
###############################################################################




###############################################################################
#
#
# Add a snipet that asks the user to check if all information looks correct and
# if so, the user needs to press a key to either continue or to break the analysis
# and correct what's wrong

echo sample sheet is found here $sample_sheet_full_location
echo flowcell_id was found to be $flowcell_id
echo will map later on the $species genome

# exit 1
#
###############################################################################




###############################################################################
#
# 0. Call bcl2fastq to demultiplex the data
#

# this needs to be optimized to automatically identify the optimal number
# of barcode-mismatches
bcl2fastq \
--no-lane-splitting --fastq-compression-level=9 \
--mask-short-adapter-reads 15 \
--barcode-mismatch 1 \
--output-dir $folder \
--sample-sheet $sample_sheet_full_location \
--runfolder-dir /data/remote/basecalls/${flowcell_id}

#
###############################################################################



    
###############################################################################
#
# 1. Call a python script that
#      - renames the fastq files to meaningful names
#      - reverses the Read1 files
#      - creates symbolic links for Read2
#

python $toolkit_folder/sequencing_preprocessing.py $folder

###############################################################################




###############################################################################
#
# 2. Start FastQC analysis
#

cd $folder

if [ ! -d fastqc ]; then
  mkdir fastqc
fi

find . -name "*.fastq.gz"  | xargs /data/rajewsky/shared_bins/FastQC-0.11.2/fastqc -t 20 -o ./fastqc/

###############################################################################




###############################################################################
#
# 3. Start the sequencing analysis pipeline
#

for i in $(find . -print | grep -i 'reversed_1.fastq.gz' | grep -v 'Undetermined');
do
    ((j=j%4)); ((j++==0)) && wait
    cd $folder;
    prefix=$(echo $i | sed 's:.*/::');
    sample_prefix=${prefix/_1.fastq.gz/};
    sample_folder=${i/$prefix/};
    cd $sample_folder;
    $toolkit_folder/sequencing_analysis.sh $species $sample_prefix &
done

wait

#
###############################################################################




###############################################################################
#
# 4. Start the pythonic script for generating the QC sheet

# first create the qc_parameters files
python $toolkit_folder/qc_sequencing_create_parameters_from_sample_sheet.py $sample_sheet_full_location

cd $folder

for i in $(find . -print | grep -i 'reversed_1.fastq.gz' | grep -v 'Undetermined');
do
    ((j=j%4)); ((j++==0)) && wait
    cd $folder;
    prefix=$(echo $i | sed 's:.*/::');
    sample_prefix=${prefix/_1.fastq.gz/};
    sample_folder=${i/$prefix/};
    cd $sample_folder$sample_prefix;
    myvar="$PWD"
    cd $toolkit_folder;
    python qc_sequencing_create_sheet.py $myvar &
done

cd $folder

#
###############################################################################
