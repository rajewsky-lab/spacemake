merged_dir = config['root_dir'] + '/projects/merged_{merged_project}/processed_data/merged_{merged_sample}/illumina/complete_data'

sample_tagged_bam = merged_dir + '/{sample}_tagged.bam'
merged_bam = merged_dir + '/merged.bam'
merged_readcounts = merged_dir + '/out_readcounts.txt.gz'
merged_qc_dir = merged_dir + qc_sheet_dir

merged_qc_sheet_parameters_file = merged_qc_dir + '/qc_sheet_parameters.yaml'

merged_reads_type_out = merged_dir + '/uniquely_mapped_reads_type.txt'
merged_top_barcodes = merged_dir + '/topBarcodes.txt'
merged_star_log_file = merged_dir + '/star_Log.final.out'

# merged dge
merged_dge_root = merged_dir + '/dge'
merged_dge_out_prefix = merged_dge_root + '/dge{dge_type}'
merged_dge_out = merged_dge_out_prefix + '.txt.gz'
merged_dge_out_summary = merged_dge_out_prefix + '_summary.txt'

#rule tag_final_bam:
#    input:
#        unpack(get_dropseq_final_bam)
#    output:
#        temporary(sample_tagged_bam)
#    threads: 4
#    shell:
#        """sambamba view -t {threads} -h {input} | \
#            awk -v suffix=.{wildcards.sample} 'BEGIN{{OFS=FS="\t"}} /^@/ {{print $0; next}} $12=$12 suffix{{print $0}}' | \
#            sambamba view -t {threads} -S /dev/stdin -o {output} -f bam"""

rule create_merged_bam:
    input:
        ancient(unpack(get_dropseq_final_bam))
    output:
        merged_bam
    threads: 4
    shell:
        "sambamba merge -t {threads} {output} {input}"

rule create_readcounts:
    input:
        merged_bam 
    output:
        merged_readcounts
    shell:
        """
        {dropseq_tools}/BamTagHistogram \
        I= {input} \
        O= {output}\
        TAG=XC
        """

rule create_merged_top_barcodes_file:
    input:
        merged_readcounts
    output:
        merged_top_barcodes
    shell:
        "set +o pipefail; zcat {input} | cut -f2 | head -100000 > {output}"

rule get_reads_type_out:
    input:
        merged_bam
    output:
        merged_reads_type_out
    shell:
        ## Script taken from sequencing_analysis.sh
        """
        samtools view {input} | \
          awk '!/GE:Z:/ && $5 == "255" && match ($0, "XF:Z:") split(substr($0, RSTART+5), a, "\t") {{print a[1]}}' | \
          awk 'BEGIN {{ split("INTRONIC INTERGENIC CODING UTR", keyword)
                      for (i in keyword) count[keyword[i]]=0
                    }}
              /INTRONIC/  {{ count["INTRONIC"]++ }}
              /INTERGENIC/  {{ count["INTERGENIC"]++ }}
              /CODING/ {{count["CODING"]++ }}
              /UTR/ {{ count["UTR"]++ }}
              END   {{
                      for (i in keyword) print keyword[i], count[keyword[i]]
                    }}' > {output}
        """

rule create_merged_dge:
    input:
        reads=merged_bam,
        top_barcodes=merged_top_barcodes
    output:
        dge=merged_dge_out,
        dge_summary=merged_dge_out_summary
    params:
        dge_root = merged_dge_root,
        dge_extra_params = lambda wildcards: get_dge_extra_params(wildcards)     
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        I= {input.reads}\
        O= {output.dge} \
        SUMMARY= {output.dge_summary} \
        CELL_BC_FILE={input.top_barcodes} \
        {params.dge_extra_params}
        """

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
        sample_params=lambda wildcards: get_qc_sheet_parameters('merged_' + wildcards.merged_sample, wildcards.umi_cutoff)
    output:
        merged_qc_sheet_parameters_file
    script:
        "qc_sequencing_create_parameters_from_sample_sheet.py"
