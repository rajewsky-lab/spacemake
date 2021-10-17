final_merged_bam = complete_data_root + final_bam_suffix + '.merged.bam'
merged_ribo_depletion_log = complete_data_root + '/ribo_depletion_log.merged.txt'
merged_star_log_file = star_prefix + 'merged.Log.final.out'

rule create_final_merged_bam:
    input:
        unpack(get_files_to_merge_snakemake(final_bam))
    output:
        final_merged_bam
    threads: 4
    shell:
        "samtools merge -n -@ {threads} -o {output} {input}"

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
        
        indices_to_save = [5, 8, 23, 10, 30]
        value_dict = {ix: 0 for ix in indices_to_save}
        indices_to_normalise = [10]

        # extract info from all logfiles, and add them up
        # we are only interested in lines 5, 8, 23, 10, 30
        # so: inp_reads, uniq_mapped_reads, avg_mapped_length, 
        # multi_mapped_reads, unmapped_too_short
        for l in logs:
            for ix in value_dict.keys():
                value_dict[ix] = value_dict[ix] + float(l[ix].split('\t')[1])

        for ix in indices_to_normalise:
            value_dict[ix] = value_dict[ix] / len(logs)

        # print to output
        with open(output[0], 'w') as fo:
            ix = 0
            for line in logs[0]:
                entry = line.split('\t') 
                if ix in value_dict.keys():
                    fo.write('%s\t%s\n' % (entry[0], value_dict[ix]))
                else:
                    fo.write('%s\t%s\n' % (entry[0], 'NA'))
                ix = ix + 1
