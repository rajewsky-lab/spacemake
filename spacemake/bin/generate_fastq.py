from spacemake.util import make_minimal_parser, setup_logging

def parse_args():
    parser = make_minimal_parser("spacemake.bin.generate_fastq")
    parser.add_argument("--CB-list", help="text file with cell barcodes to draw from")
    parser.add_argument("--UMI-n", help="length of UMI (default=16)", default=16)
    parser.add_argument("--seqs", help="text file with sequences to draw from")
    parser.add_argument("--n", help="number of times each sequence is supposed to be represented for each barcode", default=20)
    parser.add_argument("--templ1-qname", help="template string with placeholders for read1", default="@N00000:999:XXXXXXXXX:1:2101:1072:{N} 1:N:0:NNNNNNNNNN")
    parser.add_argument("--templ1-seq", help="template string with placeholders for read1", default="{CB}{UMI}")
    parser.add_argument("--templ2-qname", help="template string with placeholders for read1", default="@N00000:999:XXXXXXXXX:1:2101:1072:{N} 2:N:0:NNNNNNNNNN")
    parser.add_argument("--templ2-seq", help="template string with placeholders for read1", default="AAGCAGTGGTATCAACGCAGAGTACATGGG{seq}AAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    parser.add_argument("--out-read1", help="where to store read1", default="synthetic_R1.fastq.gz")
    parser.add_argument("--out-read2", help="where to store read2", default="synthetic_R2.fastq.gz")
    args = parser.parse_args()
    return args

def make_UMI(n, bases="ACGT"):
    import numpy as np
    umi = [bases[i] for i in np.random.randint(0, high=4, size=n)]
    return "".join(umi)

def main(args):
    import gzip
    if args.out_read1.endswith('.gz'):
        out1 = gzip.open(args.out_read1, 'wt')
    else:
        out1 = open(args.out_read1, 'w')

    if args.out_read2.endswith('.gz'):
        out2 = gzip.open(args.out_read2, 'wt')
    else:
        out2 = open(args.out_read2, 'w')

    CBs = [l.strip() for l in open(args.CB_list)]
    seqs = [l.strip() for l in open(args.seqs)]

    N = 0
    for CB in CBs:
        for seq in seqs:
            for i in range(args.n):
                UMI = make_UMI(args.UMI_n)
                qname1 = args.templ1_qname.format(**locals())
                seq1 = args.templ1_seq.format(**locals())

                qname2 = args.templ2_qname.format(**locals())
                seq2 = args.templ2_seq.format(**locals())

                qual1 = "F" * len(seq1)
                qual2 = "F" * len(seq2)

                out1.write(f"{qname1}\n{seq1}\n+\n{qual1}\n")
                out2.write(f"{qname2}\n{seq2}\n+\n{qual2}\n")
                N += 1

if __name__ == "__main__":
    args = parse_args()
    setup_logging(args, "spacemake.bin.generate_fastq")
    main(args)

    
    
    
    