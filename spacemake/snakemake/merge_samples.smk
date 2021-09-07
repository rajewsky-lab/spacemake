final_merged_bam = complete_data_root + final_bam_suffix + '.merged.bam'
merged_ribo_depletion_log = complete_data_root + '/ribo_depletion_log.merged.txt'
merged_star_log_file = star_prefix + 'merged.Log.final.out'

rule create_final_merged_bam:
    input:
        unpack(get_files_to_merge_snakemake(final_bam))
    output:
        final_merged_bam
    shell:
        "samtools merge -o {output} {input}"

rule create_merged_ribo_log:
    input:
        unpack(get_files_to_merge_snakemake(ribo_depletion_log))
    output:
        merged_ribo_depletion_log
    shell:
        "cat {input} > {output}"

rule create_merged_star_log:
    input:
        unpack(get_files_to_merge_snakemake(star_log_file))
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
