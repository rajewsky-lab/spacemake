#!/bin/bash
#  
# Master script for the computational analysis of Spatial Transcriptomics Sequencing data
#
# nikolaos.karaiskos@mdc-berlin.de
#
# version 0.1.1 
#
####################################################

####################################################
#
# Variables passed as arguments
#
# The folder where the demultiplexed data is located 
folder=$1
# The genome species that the reads will be mapped onto
species=$2
#
####################################################

####################################################
#
# 1. Call a python script that
#      - renames the fastq files to meaningful names
#      - reverses the Read1 files
#      - creates symbolic links for Read2
#

python ~/src/git/sts-sequencing/sequencing_preprocessing.py $folder

####################################################

####################################################
#
# 2. Start FastQC analysis
#

cd $folder

if [ ! -d fastqc ]; then
  mkdir fastqc
fi

find . -name "*.fastq.gz"  | xargs /data/murphy/shared_bins/FastQC-0.11.2/fastqc -t 20 -o ./fastqc/

####################################################

####################################################
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
    ~/src/git/sts-sequencing/sequencing_analysis.sh $species $sample_prefix 5000 &
done

####################################################

####################################################
#
# 4. Start the pythonic script for generating the QC sheet

# python

####################################################
