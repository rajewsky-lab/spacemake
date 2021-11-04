#!/usr/bin/env python3
__version__ = "0.9"
__author__ = ["Marvin Jens"]
__license__ = "GPL"
__email__ = ["marvin.jens@mdc-berlin.de"]

from spacemake.preprocess.fastq import parse_args, setup_logging,\
    main_combinatorial, main_dropseq

from spacemake.parallel import ExceptionLogging

def cmdline():
    with ExceptionLogging("main"):
        args = parse_args()
        NO_CALL = args.na
        setup_logging(args)

        if args.out_format == "bam" and not args.read2:
            raise ValueError("bam output format requires --read2 parameter")

        if ("bc1" in args.cell and not args.bc1_ref) or (
            "bc2" in args.cell and not args.bc2_ref
        ):
            raise ValueError(
                "bc1/2 are referenced in --cell or --cell-raw, but no reference barcodes are specified via --bc{{1,2}}-ref"
            )

        if args.bc1_ref or args.bc2_ref:
            main_combinatorial(args)
        else:
            main_dropseq(args)


if __name__ == "__main__":
    cmdline()
