#!/usr/bin/python3

#########
# about #
#########

__version__ = '0.1.1'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

###########
# imports #
###########
import os
import yaml

########
# main #
########

if __name__ == '__main__':    
    # fill the config file for the qc creation. 
    # samplesheet is parsed, so all parameters are just passed on
    params = snakemake.params['sample_params']

    with open(snakemake.output[0], 'w') as fo:

        fo.write('sample_id: ' + params['sample_id'] + '\n')
        fo.write('puck_id: ' + params['puck_id'] + '\n')
        fo.write('project_id: ' + params['project_id'] + '\n')
        fo.write('experiment: ' + params['experiment'] + '\n')
        fo.write('sequencing_date: ' + params['sequencing_date']  + '\n')
        fo.write('expected_n_beads: '+ params['input_beads'] + '\n')
        fo.write('investigator: ' + params['investigator'] + '\n')
