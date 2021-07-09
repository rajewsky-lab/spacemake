import pysam
import datetime
import argparse
import numpy as np

counted_regions = ['UTR', 'CODING']

def select_alignment(alignments):
    read_names = [aln.query_name for aln in alignments]
    if read_names.count(read_names[0]) != len(read_names):
        print(read_names)
        raise Exception(f'input alignments do not come from the same read')

    def is_exonic(aln):
        if not aln.has_tag('XF'):
            return False

        return aln.get_tag('XF') in counted_regions

    alignments_are_exonic = np.array([is_exonic(aln) for aln in alignments])

    exonic_ix = np.where(alignments_are_exonic  == True)[0]

    num_exonic = exonic_ix.shape[0]

    if num_exonic == 1:
        # if only one exonic reads from the group
        # return the exonic indices
        return alignments[exonic_ix[0]]
    else:
        return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Filter out ambiguous multi-mapper reads')

    parser.add_argument('--in-bam', help='input bam')
    parser.add_argument('--out-bam', help='output bam')

    args = parser.parse_args()
    print(args)

    bam_in = pysam.AlignmentFile(args.in_bam, "rb")

    bam_out = pysam.AlignmentFile(args.out_bam, 'wb', header= bam_in.header)
    counter = 0
    start_time = datetime.datetime.now()

    multi_mappers = []
    
    for aln in bam_in.fetch(until_eof=True):
        counter = counter + 1

        if counter % 1000000 == 0:
            finish_time = datetime.datetime.now()
            delta_seconds = (finish_time - start_time).seconds
            # restart time
            start_time = finish_time
            print(f'Processed 1 millon records in {delta_seconds} seconds, total records processed {counter}. current time: {finish_time}')

        mapped_number = aln.get_tag('NH')

        if mapped_number == 1:
            bam_out.write(aln)
        else:
            if len(multi_mappers) < (mapped_number - 1):
                # still some multimappers missing. we need to add the alignments 
                # until the last one to the list
                multi_mappers.append(aln)
            else:
                # add the last alignment
                multi_mappers.append(aln)
                # decide which, if any, to keep
                aln_to_keep = select_alignment(multi_mappers)

                if aln_to_keep is not None:
                    # set aln secondary flag to 0, so that it is flagged as primary
                    # secondary flag is at 0x100, so 8th bit (starting from 0)
                    aln_to_keep.flag = aln_to_keep.flag & ~(1<<8)
                    bam_out.write(aln_to_keep)
                    
                # reset multimapper list
                multi_mappers = []
