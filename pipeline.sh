#!/bin/bash
#  
# Master pipeline script for the computational analysis of Spatial Transcriptomics Sequencing data
#
# author : Nikos Karaiskos
# email : nikolaos.karaiskos@mdc-berlin.de
#
# version 0.1.3
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
# 3. Start the dropseq pipeline tailored for sts-seq
#

for i in $(find . -print | grep -i 'reversed_1.fastq.gz');
do
    ((j=j%4)); ((j++==0)) && wait
    cd $folder;
    prefix=$(echo $i | sed 's:.*/::');
    sample_prefix=${prefix/_1.fastq.gz/};
    sample_folder=${i/$prefix/};
    cd $sample_folder;
    $toolkit_folder/sequencing_analysis.sh $species $sample_prefix &
done

###############################################################################

###############################################################################
#
# 4. Start the pythonic script for generating the QC sheet

for i in $(find . -print | grep -i 'reversed_1.fastq.gz');
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
    echo "Done";
done

###############################################################################
