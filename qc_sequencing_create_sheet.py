#!/usr/bin/python3

#########
# about #
#########

__version__ = '0.1.11'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

###########
# imports #
###########
import time
import gzip
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fpdf
from collections import Counter
import sys
import yaml
import subprocess
import itertools
import math
from datetime import datetime

#############
# snakemake #
#############
# 
# This file will be called by snakemake. 
# The variable snakemake holds all the input data which is needed to run functions here
#
#############
# functions #
#############

def is_gzip_file(filename):
    """Test if a file is gzip compressed or not."""
    try:
        # This will raise OSError for uncompressed files & has no side
        # effects for compressed files:
        gzip.GzipFile(filename).peek(1)
        return True
    except OSError:
        return False

def open_file(filename):
    """Open a file (gzip compressed or not)."""
    if (is_gzip_file(filename)):
        return gzip.open(filename,'rt')
    else:
        return open(filename,'rt')

def find_filename(folder, endswith=''):
    """Find a filename if its ending is known. Assumes that the ending is unique.
    folder   -- The folder containing the file
    endswith -- The string that the filename ends with, including extension."""
    filename = [os.path.join(root, f) for root, _, files in os.walk(folder)
                        for f in files
                        if f.endswith(endswith)][0]
    return filename

def compute_shannon_entropy(barcode):
    prob, length = Counter(barcode), float(len(barcode))
    return -sum( count/length * math.log(count/length, 2) for count in prob.values())


def compress_string(barcode):
    return ''.join(
            letter + str(len(list(group)))
            for letter, group in itertools.groupby(barcode))

def load_read_statistics():
    read_statistics = dict()

    with open(snakemake.input.star_log, 'r') as fi:
        idx = 0
        for line in fi.readlines():
            entry = line.strip('\n').split('\t')
            if idx == 5:
                read_statistics['input reads'] = int(entry[1])
            if idx == 8:
                read_statistics['uniquely mapped'] = int(entry[1])
            idx += 1
 
    read_types = pd.read_csv(snakemake.input.reads_type_out, names=['name', 'num'], sep=' ') 

    read_statistic_keys = ['intronic', 'intergenic', 'amb']

    for key in read_statistic_keys:
        read_statistics[key] = 0

        # get the read statistic count as a list. list has 1 element if staistic exist, 0 if it doesnt
        count = read_types.num[read_types.name == key.upper()].to_list()

        if len(count) == 1:
            read_statistics[key] = count[0]

    return read_statistics

def load_bead_statistics(folder):
    """Read basic stastistics concerning the beads
    folder -- The folder containing the sequencing data after analyzing it
              with the standard Dropseq pipeline."""

    # the dictionary holding the values
    bead_statistics = dict()
    
    # read readcounts for all barcodes seen in the data
    readcounts = pd.read_csv(snakemake.input.read_counts, sep='\t',
        skiprows=1, names=['reads', 'barcode'])
    bead_statistics['total # of barcodes'] = readcounts.shape[0]

    # select # barcodes for cumulative fraction
    barcode_limit = min(100000, bead_statistics['total # of barcodes'])
    rc_cumsum = readcounts['reads'][:barcode_limit].cumsum()
    rc_cumsum /= max(rc_cumsum)
    plt.plot(np.arange(barcode_limit), rc_cumsum)
    plt.xlabel('bead barcodes sorted by number of reads', fontsize=18)
    plt.ylabel('cumulative fraction of reads', fontsize=18)
    plt.savefig(folder+'cumulative_fraction.png')
    plt.tight_layout()
    plt.close()

    # plot read number distribution for the first 100,000 beads
    bead_reads = np.array(readcounts['reads'][:100000])
    plt.hist(bead_reads, bins=1000)
    plt.xlabel('# reads / bead', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(folder+'bead_reads_distribution.png')
    plt.close()

    # calculate Shannon entropies for the barcodes
    barcode_entropies = np.round(np.array([compute_shannon_entropy(bc) for 
        bc in readcounts['barcode'][:barcode_limit]]), 2)
    bead_statistics['barcode_entropies'] = barcode_entropies
    plt.hist(barcode_entropies, bins=100)
    plt.xlabel('Shannon entropy', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.tight_layout()
    plt.savefig(folder+'barcode_entropies.png')
    plt.close()

    # calculate string compression for the barcodes
    barcode_string_compr = np.array([len(compress_string(bc)) for 
        bc in readcounts['barcode'][:barcode_limit]])
    bead_statistics['barcode_string_compression'] = barcode_string_compr
    plt.hist(barcode_string_compr, bins=20)
    plt.xlabel('string compression', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.tight_layout()
    plt.savefig(folder+'barcode_string_compression.png')
    plt.close()

    return bead_statistics

def load_downstream_statistics(folder, umi_cutoff):
    """Read the load_downstream_statistics.
    folder    -- The folder containing the sequencing data after analyzing it
              with the standard Dropseq pipeline.
    umi_cutoff -- The minimum number of UMIs to keep a bead."""

    # the dictionary holding the values
    downstream_statistics = dict()

    # read the summary table of dge_all (containing intronic and exonic reads)
    # change this back to snakemake output before merging
    downstream_stats = pd.read_csv(snakemake.input.dge_all_summary,
            # skip the first 5 rows as they contain comments
            skiprows=6,
            # rename the column, make column=0 the index
            sep='\t', index_col=0).rename(columns={'NUM_GENIC_READS': 'reads', 'NUM_TRANSCRIPTS':'umis', 'NUM_GENES':'genes'})

    # find bead number distribution as a function of the umi threshold
    bead_dist = np.array([np.sum(np.array(downstream_stats['umis']) >= x) for x in range(100)])
    # and plot this for the qc_sheet
    plt.plot(np.arange(100), bead_dist)
    plt.xlabel('UMI threshold', fontsize=18)
    plt.ylabel('# beads', fontsize=18)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(folder+'bead_distribution_umi_threshold.png')
    plt.close()

    # we decrease the treshold if the output would be empty otherwise
    while(sum(downstream_stats['umis'] >= umi_cutoff) == 0):
        umi_cutoff = umi_cutoff - 10

    downstream_statistics['minimum umis per bead'] = umi_cutoff

    # filter by umi_cutoff given in the qc_sequencing_parameters.yaml
    downstream_stats = downstream_stats[downstream_stats['umis'] >= umi_cutoff]
    
    # find beads which have the minimum number of UMIs
    beads = downstream_stats.index.str.split('.').str[0].to_list()
    downstream_statistics['beads'] = len(beads)

    # compute total reads per bead and plot histogram
    reads_per_bead = downstream_stats['reads']
    downstream_statistics['reads per bead'] = int(round(reads_per_bead.median()))
    ax = reads_per_bead.hist(bins=100, ylabelsize=18, xlabelsize=18)
    plt.xlabel('reads per bead', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(folder+'hist_reads_per_bead.png')
    plt.close()

    # compute total genes per bead and plot histogram
    genes_per_bead = downstream_stats['genes']
    downstream_statistics['genes per bead'] = int(round(genes_per_bead.median()))
    ax = genes_per_bead.hist(bins=100, ylabelsize=18, xlabelsize=18)
    plt.xlabel('genes per bead', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(folder+'hist_genes_per_bead.png')
    plt.close()

    # compute total umis per bead and plot histogram
    umis_per_bead = downstream_stats['umis']
    downstream_statistics['umis per bead'] = round(int(umis_per_bead.median()))
    ax = umis_per_bead.hist(bins=100, ylabelsize=18, xlabelsize=18)
    plt.xlabel('umis per bead', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(folder+'hist_umis_per_bead.png')
    plt.close()

    # split barcodes to individual bases
    bases = np.array([x for y in range(len(beads)) for x in beads[y]])
    # bast to 2D array
    beads = bases.reshape(len(beads), len(beads[0]))
    beads = beads.T
    # count the nucleotide frequencies
    dict_list = [dict(Counter(beads[x])) for x in range(beads.shape[0])]
    nt_composition = pd.DataFrame(dict_list)
    # make a plot and save it on the disk
    nt_composition.plot.bar()
    plt.legend(loc='lower right')
    plt.xlabel('barcode position', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.tight_layout()
    plt.savefig(folder+'nucleotide_composition.png')

    return downstream_statistics

def create_qc_sheet(folder):
    # Uncomment the following 2 lines before merging
    with open(snakemake.input.parameters_file) as f:
        parameters = yaml.load(f, Loader=yaml.FullLoader)

    # Comment these lines before merging (keep them for future enhancements)
    # parameters = dict()
    # parameters['umi_cutoff'] = 100
    # parameters['project_id'] = 'some_id'
    # parameters['sample_id'] = 'some_id'
    # parameters['puck_id'] = 'some_id'
    # parameters['experiment'] = 'some_experiment'
    # parameters['sequencing_date'] = 'some_date'
    # parameters['investigator'] = 'some_name'
    # parameters['input_beads'] = '100000'    

    read_statistics = load_read_statistics()
    bead_statistics = load_bead_statistics(folder)
    downstream_statistics = load_downstream_statistics(folder, umi_cutoff=parameters['umi_cutoff'])

    input_reads = read_statistics['input reads']
    uniq_mapped = read_statistics['uniquely mapped']
    intronic = read_statistics['intronic']
    intergenic = read_statistics['intergenic']
    ambiguous = read_statistics['amb']

    pdf = fpdf.FPDF()
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('Arial', '', 11)
    pdf.cell(90, 10, " ", 0, 1, 'C')
    pdf.cell(90, 10, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.cell(100, 5, "project_id: %s" % (parameters['project_id']), 0, 2, 'L')
    pdf.cell(100, 5, "sample_id: %s" % (parameters['sample_id']), 0, 2, 'L')
    pdf.cell(100, 5, "puck_id: %s" % (parameters['puck_id']), 0, 2, 'L')
    # sometimes we have more experiment descriptions. this is when we have merged samples.
    # in this case we need to print each experiment desc, one below the other
    # first just print one cell titled experiment
    pdf.cell(25, 5, "experiment: ", 0, 0, 'L')
    # then right of it print one cell per experiment description. for merged they are split by ','
    for experiment in parameters['experiment'].split(','):
        pdf.cell(75, 5, experiment, 0, 2, 'L')
    # create an empty cell to reposition
    pdf.cell(75, 0, ' ', 0, 1, 'L')
    pdf.cell(10)
    pdf.cell(100, 5, "sequencing_date: %s" % (', '.join(parameters['sequencing_date'].split(','))), 0, 2, 'L')
    pdf.cell(100, 5, "investigator(s): %s" % (', '.join(parameters['investigator'].split(','))), 0, 2, 'L')
    pdf.cell(90, 5, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.cell(100, 5, "sequencing QC v."+str(__version__)+ 
        ", nikolaos.karaiskos@mdc-berlin.de, tamasryszard.sztanka-toth@mdc-berlin.de", 0, 2, 'L')
    pdf.cell(100, 5, "QC generated on %s" % (datetime.now().strftime('%d/%m/%Y %H:%M')), 0, 2, 'L')
    pdf.cell(90, 8, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.set_font('Arial', '', 10)
    pdf.cell(35, 8, 'input reads', 1, 0, 'C')
    pdf.cell(35, 8, 'uniquely mapped', 1, 0, 'C')
    pdf.cell(35, 8, 'intergenic', 1, 0, 'C')
    pdf.cell(35, 8, 'intronic', 1, 0, 'C')
    pdf.cell(35, 8, 'ambiguous', 1, 1, 'C')
    pdf.cell(10)
    pdf.cell(35, 8, format(input_reads, ','), 1, 0, 'C')
    pdf.cell(35, 8, format(uniq_mapped, ',')
        + ' (' + str(round(100*uniq_mapped/input_reads, 1)) + '%)' , 1, 0, 'C')
    pdf.cell(35, 8, format(intergenic, ',') 
        + ' (' + str(round(100*intergenic/uniq_mapped, 1)) + '%)', 1, 0, 'C')
    pdf.cell(35, 8, format(intronic, ',') 
        + ' (' + str(round(100*intronic/uniq_mapped, 1)) + '%)', 1, 0, 'C')
    pdf.cell(35, 8, format(ambiguous, ',') 
        + ' (' + str(round(100*ambiguous/uniq_mapped, 1)) + '%)', 1, 1, 'C')
    pdf.cell(90, 5, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.cell(35, 8, 'input # beads', 1, 0, 'C')
    pdf.cell(35, 8, 'total # beads', 1, 1, 'C')
    pdf.cell(10)
    pdf.cell(35, 8, str(parameters['input_beads']), 1, 0, 'C')
    pdf.cell(35, 8, format(bead_statistics['total # of barcodes'], ','), 1, 1, 'C')
    pdf.cell(90, 5, " ", 0, 2, 'C')
    pdf.cell(10)
    pdf.set_font('Arial', '', 11)
    pdf.image(folder+'cumulative_fraction.png', x=None, y=None, w=75, h=50, type='', link='')
    pdf.set_xy(pdf.get_x()+85, pdf.get_y()-46)
    pdf.image(folder+'bead_reads_distribution.png', x=100, y=118, w=75, h=50, type='', link='')
    pdf.image(folder+'bead_distribution_umi_threshold.png', x=20, y=180, w=75, h=50, type='', link='')
    
    # 2nd page
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.cell(90, 15, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.cell(20, 8, 'beads', 1, 0, 'C')
    pdf.cell(20, 8, 'reads', 1, 0, 'C')
    pdf.cell(20, 8, 'genes', 1, 0, 'C')
    pdf.cell(20, 8, 'umis', 1, 0, 'C')
    pdf.cell(25, 8, 'min umi cutoff', 1, 1, 'C')
    pdf.cell(10)
    pdf.cell(20, 8, format(downstream_statistics['beads'], ','), 1, 0, 'C')
    pdf.cell(20, 8, format(downstream_statistics['reads per bead'], ','), 1, 0, 'C')
    pdf.cell(20, 8, format(downstream_statistics['genes per bead'], ','), 1, 0, 'C')
    pdf.cell(20, 8, format(downstream_statistics['umis per bead'], ','), 1, 0, 'C')
    pdf.cell(25, 8, format(downstream_statistics['minimum umis per bead'], ','), 1, 1, 'C')
    pdf.cell(90, 8, " ", 0, 2, 'C')
    pdf.image(folder+'hist_reads_per_bead.png', x=20, y=50, w=75, h=50, type='', link='')
    pdf.image(folder+'hist_genes_per_bead.png', x=100, y=50, w=75, h=50, type='', link='')
    pdf.image(folder+'hist_umis_per_bead.png', x=20, y=110, w=75, h=50, type='', link='')
    pdf.image(folder+'barcode_entropies.png', x=100, y=110, w=75, h=50, type='', link='')
    pdf.image(folder+'barcode_string_compression.png', x=20, y=170, w=75, h=50, type='', link='')
    pdf.image(folder+'nucleotide_composition.png', x=100, y=170, w=75, h=50, type='', link='')

    sample_folder = folder.strip('/$').split('/')
    sample_folder = sample_folder[-1]
    pdf.output(snakemake.output[0], 'F') # uncomment before merging


########
# main #
########

if __name__ == '__main__':
    
    start_time = time.time()

    # provide the folder as an argument
    qc_sheet_folder = os.path.dirname(snakemake.output[0]) + '/'
    print ('starting analysis for sample in folder', qc_sheet_folder, '... ')

    # subprocess.call('mkdir -p ' + qc_sheet_folder, shell=True) # uncomment before merging
    
    create_qc_sheet(qc_sheet_folder)

    print ('took', round(time.time() - start_time, 2), 'seconds')
