#!/usr/bin/python3

#########
# about #
#########

__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de']

###########
# imports #
###########
import os
import yaml
import sys

#############
# Variables #
#############
toolkit_folder = '~/src/git/sts-sequencing/'

#############
# functions #
#############

def create_qc_sequencing_parameters(sample_sheet):
    return

########
# main #
########

if __name__ == '__main__':    

    # provide the sample_sheet location as an argument
    # it assumes that the data is already demultiplexed and therefore all
    # relevant sample directories have been created
    sample_sheet = sys.argv[1]
    sample_sheet_folder = '/'.join(sys.argv[1].split('/')[:-1])

    with open(sample_sheet, 'r') as fi:
        sample_id_seen = False
        for line in fi:
            line = line.strip('\n')
            if 'Investigator' in line:
                investigator = line.split(',')[1]
            if 'Date' in line:
                sequencing_date = line.split(',')[1]
            if 'Sample_ID' in line:
                sample_id_seen = True
                continue
            if sample_id_seen:
                sample_info = line.split(',')
                if len(sample_info) == 1:
                    break

                # print (sample_info[6] + '/' + sample_info[0] + 
                #     '/' + sample_info[1])
                # print ('file extension for the qc_sheet in the end',
                #     sample_info[0] + '_' + sample_info[1])

                with open(sample_sheet_folder + '/demultiplexed_data/' + 
                          sample_info[6] + '/' + sample_info[0] + '/' +
                          sample_info[1] + '_reversed/output_qc_sheet/' +
                          'qc_sequencing_parameters_' + sample_info[0] + 
                          '_' + sample_info[1] + '.yaml', 'w') as fo:

                    fo.write('sample_id: ' + sample_info[0] + '\n')
                    fo.write('puck_id: ' + sample_info[1] + '\n')
                    fo.write('experiment: ' + sample_info[7] + '\n')
                    fo.write('date: ' + sequencing_date + '\n')
                    fo.write('input_beads: 60k-100k' + '\n')
                    fo.write('threshold: 100' + '\n')








    












