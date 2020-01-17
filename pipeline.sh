#!/bin/bash
#  
# Master pipeline script for the computational analysis of Spatial Transcriptomics Sequencing data
#
# author : Nikos Karaiskos
# email : nikolaos.karaiskos@mdc-berlin.de
#
# version 0.1.5
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
# The folder where the demultiplexed data is located 
folder=$1

# The genome species that the reads will be mapped onto
species=$2

# The sample sheet. right now it's the 3rd parameter but eventually
# will replace $1 so that the demultiplexing starts there
sample_sheet=$3

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

###############################################################################

###############################################################################
#
# 4. Start the pythonic script for generating the QC sheet

# first create the qc_parameters files
python $toolkit_folder/qc_sequencing_create_parameters_from_sample_sheet.py $sample_sheet

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

###############################################################################
