#!/usr/bin/env python3
from spacemake.contrib import __version__, __license__, __author__, __email__
from spacemake.preprocess.fastq import (
    parse_args,
    main_combinatorial,
    main_dropseq,
    # main_bulk
)
from spacemake.util import setup_logging
from spacemake.parallel import ExceptionLogging


def cmdline():
    args = parse_args()
    NO_CALL = args.na
    setup_logging(args, name="spacemake.preprocess")

    with ExceptionLogging("spacemake.preprocess.main") as el:
        if args.out_format == "bam" and not args.read2:
            raise ValueError("bam output format requires --read2 parameter")

        if ("bc1" in args.cell and not args.bc1_ref) or (
            "bc2" in args.cell and not args.bc2_ref
        ):
            raise ValueError(
                "bc1/2 are referenced in --cell or --cell-raw, but no reference barcodes are specified via --bc{{1,2}}-ref"
            )

        if args.bc1_ref or args.bc2_ref:
            res = main_combinatorial(args)
        else:
            res = main_dropseq(args)

    if el.exception:
        return -1
    else:
        el.logger.info(f"exit code={res}")
        return res


if __name__ == "__main__":
    import sys

    ret_code = cmdline()
    sys.exit(ret_code)
