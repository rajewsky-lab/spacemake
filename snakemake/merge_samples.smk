merged_dir = config['root_dir'] + '/projects/merged_{merged_project}/processed_data/merged_{merged_sample}/illumina/complete_data'

sample_tagged_bam = merged_dir + '/{sample}_tagged.bam'
merged_bam = merged_dir + '/merged.bam'
merged_final_bam = merged_dir + '/final.bam'

merged_qc_sheet_parameters_file = merged_dir + '/qc_sheet_parameters.yaml'

merged_star_log_file = merged_dir + '/star_Log.final.out'
merged_ribo_depletion_log = merged_dir + '/ribo_depletion_log.txt'

rule create_merged_bam:
    input:
        ancient(unpack(get_dropseq_final_bam))
    output:
        merged_bam
    threads: 4
    shell:
        "sambamba merge -t {threads} {output} {input}"

rule create_final_bam:
    input:
        merged_bam
    output:
        merged_final_bam
    shell:
       'ln -sr {input} {output}' 

rule create_merged_ribo_depletion_log:
    input:
        unpack(get_merged_ribo_depletion_log_inputs)
    output:
        merged_ribo_depletion_log
    shell:
        "cat {input} > {output}"

rule create_merged_star_log:
    input:
        unpack(get_merged_star_log_inputs)
    output:
        merged_star_log_file
    run:
        logs = []
        for f in input:
            with open(f, 'r') as fi:
                logs = logs + [fi.read().splitlines()]

        inp_reads = 0
        uniquely_mapped = 0

        # extract info from all logfiles, and add them up
        for l in logs:
            inp_reads = inp_reads + int(l[5].split('\t')[1])
            uniquely_mapped = uniquely_mapped + int(l[8].split('\t')[1])

        # print to output
        with open(output[0], 'w') as fo:
            idx = 0
            for line in logs[0]:
                entry = line.split('\t') 
                if idx == 5:
                    fo.write('%s\t%s\n' % (entry[0], inp_reads))
                elif idx == 8:
                    fo.write('%s\t%s\n' % (entry[0], uniquely_mapped))
                else:
                    fo.write('%s\t%s\n' % (entry[0], 'NA'))
                idx = idx + 1

rule create_merged_qc_parameters:
    params:
        sample_params=lambda wildcards: get_qc_sheet_parameters('merged_' + wildcards.project, 'merged_' + wildcards.sample)
    output:
        merged_qc_sheet_parameters_file
    script:
        "analysis/qc_sequencing_create_parameters_from_sample_sheet.py"
