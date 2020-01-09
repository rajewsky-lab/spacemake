#!/usr/bin/python3

#########
# about #
#########

__version__ = '0.1.7'
__author__ = ['Nikos Karaiskos']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de']

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
import collections
import sys
import yaml
import subprocess

#############
# Variables #
#############
toolkit_folder = '~/src/git/sts-sequencing/'

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

def load_read_statistics(folder):
    """Read the basic and mapped read stastistics.
    folder -- The folder containing the sequencing data after analyzing it
              with the standard Dropseq pipeline."""

    read_statistics = dict()

    with open(folder + 'star_Log.final.out', 'r') as fi:
        idx = 0
        for line in fi.readlines():
            entry = line.strip('\n').split('\t')
            if idx == 5:
                read_statistics['input reads'] = int(entry[1])
            if idx == 8:
                read_statistics['uniquely mapped'] = int(entry[1])
            idx += 1

    with open(folder + 'uniquely_mapped_reads_type.txt', 'r') as fi:
        idx = 0
        for line in fi.readlines():
            entry = line.strip('\n').split(' ')
            if idx == 1:
                read_statistics['intronic'] = int(entry[1])
            if idx == 2:
                read_statistics['intergenic'] = int(entry[1])
            idx += 1

    return read_statistics

def load_bead_statistics(folder):
    """Read basic stastistics concerning the beads
    folder -- The folder containing the sequencing data after analyzing it
              with the standard Dropseq pipeline."""

    # the dictionary holding the values
    bead_statistics = dict()
    
    # read readcounts for all barcodes seen in the data
    readcounts = pd.read_csv(folder + 'out_readcounts.txt.gz', sep='\t', 
        skiprows=1, names=['reads', 'barcode'])
    bead_statistics['total # of barcodes'] = readcounts.shape[0]

    # select # barcodes for cumulative fraction
    barcode_limit = 25000
    rc_cumsum = readcounts['reads'][:barcode_limit].cumsum()
    rc_cumsum /= max(rc_cumsum)
    plt.plot(np.arange(barcode_limit), rc_cumsum)
    plt.xlabel('bead barcodes sorted by number of reads', fontsize=18)
    plt.ylabel('cumulative fraction of reads', fontsize=18)
    plt.savefig(folder+'output_qc_sheet/cumulative_fraction.png')
    plt.tight_layout()
    plt.close()

    # read the synthesis errors summary from the dropseq toolkit
    filename = find_filename(folder, 'synthesis_stats_summary.txt')
    with open(filename, 'r') as fi:
        idx = 0
        for line in fi.readlines():
            if idx == 3:
                entry = line.strip('\n').split('\t')
                break
            idx += 1
    pct = int(entry[1]) / int(entry[0])
    pct = round(pct * 100, 2)
    bead_statistics['beads without synthesis errors'] = pct

    # read the substitution errors file from the dropseq toolkit
    filename = find_filename(folder, 'substitution_errors_report.txt')
    with open(filename, 'r') as fi:
        idx = 0
        for line in fi.readlines():
            entry = line.strip('\n').split('=')
            if idx == 5:
                total_barcodes_tested = int(entry[1])
            if idx == 6:
                barcodes_collapsed = int(entry[1])
            idx += 1
    bead_statistics['barcodes collapsed'] = barcodes_collapsed

    return bead_statistics

def load_downstream_statistics(folder, threshold):
    """Read the load_downstream_statistics.
    folder    -- The folder containing the sequencing data after analyzing it
              with the standard Dropseq pipeline.
    threshold -- The minimum number of UMIs to keep a bead."""

    # the dictionary holding the values
    downstream_statistics = dict()
    downstream_statistics['minimum umis per bead'] = threshold

    # read the Digital Gene Expression matrix
    start_time = time.time()
    print ('Loading the DGE... ', end='', flush=True)

    # call the R script that calculates the downstream statistics
    subprocess.call("Rscript " + toolkit_folder + 
                    "qc_sequencing_generate_downstream_statistics.R " + folder + ' '
                    + str(threshold), shell=True)
    # read the result of the Rscript here through pandas
    downstream_stats_R = pd.read_csv(folder+'output_qc_sheet/downstream_statistics.csv',
                                     index_col=0)
    
    print ('[', round(time.time()-start_time, 2), 'seconds ]')

    # find beads which have the minimum number of UMIs
    beads = downstream_stats_R.index.tolist()
    downstream_statistics['beads'] = len(beads)

    # compute total reads per bead and plot histogram
    reads_per_bead = downstream_stats_R['reads']
    downstream_statistics['reads per bead'] = int(round(reads_per_bead.median()))
    ax = reads_per_bead.hist(bins=100, ylabelsize=18, xlabelsize=18)
    plt.xlabel('reads per bead', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(folder+'output_qc_sheet/hist_reads_per_bead.png')
    plt.close()

    # compute total genes per bead and plot histogram
    genes_per_bead = downstream_stats_R['genes']
    downstream_statistics['genes per bead'] = int(round(genes_per_bead.median()))
    ax = genes_per_bead.hist(bins=100, ylabelsize=18, xlabelsize=18)
    plt.xlabel('genes per bead', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(folder+'output_qc_sheet/hist_genes_per_bead.png')
    plt.close()

    # compute total umis per bead and plot histogram
    umis_per_bead = downstream_stats_R['umis']
    downstream_statistics['umis per bead'] = round(int(umis_per_bead.median()))
    ax = umis_per_bead.hist(bins=100, ylabelsize=18, xlabelsize=18)
    plt.xlabel('umis per bead', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.xscale('log')
    plt.tight_layout()
    plt.savefig(folder+'output_qc_sheet/hist_umis_per_bead.png')
    plt.close()

    # split barcodes to individual bases
    bases = np.array([x for y in range(len(beads)) for x in beads[y]])
    # bast to 2D array
    beads = bases.reshape(len(beads), len(beads[0]))
    beads = beads.T
    # count the nucleotide frequencies
    dict_list = [dict(collections.Counter(beads[x])) for x in range(beads.shape[0])]
    nt_composition = pd.DataFrame(dict_list)
    # make a plot and save it on the disk
    nt_composition.plot.bar()
    plt.legend(loc='lower right')
    plt.xlabel('barcode position', fontsize=18)
    plt.ylabel('count', fontsize=18)
    plt.tight_layout()
    plt.savefig(folder+'output_qc_sheet/nucleotide_composition.png')

    return downstream_statistics

def create_qc_sheet(folder):
    with open(folder+'qc_sequencing_parameters.yaml') as f:
        parameters = yaml.load(f, Loader=yaml.FullLoader)

    read_statistics = load_read_statistics(folder)
    bead_statistics = load_bead_statistics(folder)
    downstream_statistics = load_downstream_statistics(folder, threshold=parameters['threshold'])

    input_reads = read_statistics['input reads']
    uniq_mapped = read_statistics['uniquely mapped']
    intronic = read_statistics['intronic']
    intergenic = read_statistics['intergenic']

    pdf = fpdf.FPDF()
    pdf.add_page()
    pdf.set_xy(0, 0)
    pdf.set_font('Arial', '', 11)
    pdf.cell(90, 10, " ", 0, 1, 'C')
    pdf.cell(90, 10, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.cell(100, 8, ", ".join([parameters['puck_id'], 
                                parameters['experiment'], 
                                str(parameters['date'])]), 0, 1, 'L')
    pdf.cell(10)
    pdf.cell(100, 8, "sequencing QC v."+str(__version__)+ 
        ", nikolaos.karaiskos@mdc-berlin.de", 0, 1, 'L')
    pdf.cell(90, 8, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.cell(30, 8, 'input reads', 1, 0, 'C')
    pdf.cell(40, 8, 'uniquely mapped', 1, 0, 'C')
    pdf.cell(40, 8, 'intergenic', 1, 0, 'C')
    pdf.cell(40, 8, 'intronic', 1, 1, 'C')
    pdf.cell(10)
    pdf.cell(30, 8, format(input_reads, ','), 1, 0, 'C')
    pdf.cell(40, 8, format(uniq_mapped, ',')
        + ' (' + str(round(100*uniq_mapped/input_reads, 1)) + '%)' , 1, 0, 'C')
    pdf.cell(40, 8, format(intergenic, ',') 
        + ' (' + str(round(100*intergenic/uniq_mapped, 1)) + '%)', 1, 0, 'C')
    pdf.cell(40, 8, format(intronic, ',') 
        + ' (' + str(round(100*intronic/uniq_mapped, 1)) + '%)', 1, 1, 'C')
    pdf.cell(90, 5, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.cell(30, 8, 'input # beads', 1, 0, 'C')
    pdf.cell(30, 8, 'total # beads', 1, 0, 'C')
    pdf.cell(60, 8, 'beads without synth errors', 1, 0, 'C')
    pdf.cell(30, 8, 'bc collapsed', 1, 1, 'C')
    pdf.cell(10)
    pdf.cell(30, 8, str(parameters['input_beads']), 1, 0, 'C')
    pdf.cell(30, 8, format(bead_statistics['total # of barcodes'], ','), 1, 0, 'C')
    pdf.cell(60, 8, str(bead_statistics['beads without synthesis errors'])
        +  '%', 1, 0, 'C')
    pdf.cell(30, 8, format(bead_statistics['barcodes collapsed'], ','), 1, 1, 'C')
    pdf.cell(90, 5, " ", 0, 2, 'C')
    pdf.cell(10)
    pdf.image(folder+'output_qc_sheet/cumulative_fraction.png', x=None, y=None, w=75, h=50, type='', link='')
    pdf.image(folder+'output_qc_sheet/nucleotide_composition.png', x=100, y=89, w=75, h=50, type='', link='')
    pdf.cell(90, 5, " ", 0, 1, 'C')
    pdf.cell(10)
    pdf.cell(20, 8, 'beads', 1, 0, 'C')
    pdf.cell(20, 8, 'reads', 1, 0, 'C')
    pdf.cell(20, 8, 'genes', 1, 0, 'C')
    pdf.cell(20, 8, 'umis', 1, 0, 'C')
    pdf.cell(25, 8, 'threshold', 1, 1, 'C')
    pdf.cell(10)
    pdf.cell(20, 8, format(downstream_statistics['beads'], ','), 1, 0, 'C')
    pdf.cell(20, 8, format(downstream_statistics['reads per bead'], ','), 1, 0, 'C')
    pdf.cell(20, 8, format(downstream_statistics['genes per bead'], ','), 1, 0, 'C')
    pdf.cell(20, 8, format(downstream_statistics['umis per bead'], ','), 1, 0, 'C')
    pdf.cell(25, 8, format(downstream_statistics['minimum umis per bead'], ','), 1, 1, 'C')
    pdf.cell(90, 5, " ", 0, 2, 'C')
    pdf.set_font('arial', '', 12)
    pdf.image(folder+'output_qc_sheet/hist_reads_per_bead.png', x=None, y=None, w=75, h=50, type='', link='')
    pdf.image(folder+'output_qc_sheet/hist_genes_per_bead.png', x=100, y=162, w=75, h=50, type='', link='')
    pdf.image(folder+'output_qc_sheet/hist_umis_per_bead.png', x=None, y=None, w=75, h=50, type='', link='')
    pdf.output(folder+'output_qc_sheet/qc_sheet.pdf', 'F')


########
# main #
########

if __name__ == '__main__':
    
    start_time = time.time()

    # provide the folder as an argument
    folder = '/' + sys.argv[1].strip('/') + '/'
    print ('starting analysis for sample in folder', folder, '... ')

    try:
        os.mkdir(folder + 'output_qc_sheet')
    except OSError:
        print ("Skipped creating new output folder, already exists")
    
    create_qc_sheet(folder)

    print ('took', round(time.time() - start_time, 2), 'seconds')
