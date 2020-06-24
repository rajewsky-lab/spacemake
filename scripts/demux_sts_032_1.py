# trying to rescue reads from undetermined demultiplex fastq.gz.
# correct barcode is TAAGGCGA, for the ercc (sts_032_1) sample
# if 2nd and3rd nt is AA, and nt6 is C or nt8 is A, then we can rescue it

import gzip

line_c = 0
barcode = 'TAAGGCGA'

rescue_read = False

with gzip.open('/data/rajewsky/projects/slide_seq/demultiplex_data/200617_sts_032_sts_K_003_3_4_sample_sheet/Undetermined_S0_R1_001.fastq.gz', 'rt') as r1,\
    gzip.open('/data/rajewsky/projects/slide_seq/demultiplex_data/200617_sts_032_sts_K_003_3_4_sample_sheet/Undetermined_S0_R2_001.fastq.gz', 'rt') as r2,\
    gzip.open('out/sts_032_1_extra_R1.fastq.gz', 'wt') as r1_out, gzip.open('out/sts_032_1_extra_R2.fastq.gz', 'wt') as r2_out:
        for line1, line2 in zip(r1, r2):
            if line_c % 4 == 0:
                # get the real demux barcode. also remove the \n from the end
                demux_barcode = line1.split(':')[-1][:8]

                if demux_barcode[1] == barcode[1] and demux_barcode[2] == barcode[2] and (demux_barcode[5] == barcode[5] or demux_barcode[7] == barcode[7]):
                    rescue_read = True
                else:
                    rescue_read = False

            if rescue_read:
                r1_out.write(line1)
                r2_out.write(line2)

            line_c = line_c + 1
