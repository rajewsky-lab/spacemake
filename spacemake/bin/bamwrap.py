import mrfifo as mf


def parse_args():
    import spacemake.util as util

    parser = util.make_minimal_parser(
        "bamwrap",
        description="pipe a BAM file into a tool (aligner?) that expects FASTA/Q and outputs SAM. Merge BAM tags after the fact.",
    )
    parser.add_argument(
        "--bam-in",
        default="/dev/stdin",
        help="bam input (default=/dev/stdin)",
    )
    parser.add_argument(
        "--format",
        default="fastq",
        choices=["fastq", "fasta"],
        help="what to extract from the BAM. fastq or fasta (default=fastq)",
    )
    parser.add_argument(
        "--cmd",
        help="quote-encapsulated commandline you wish to run, which must use {input} and {ouput}, and treat these as files",
        default="cat {input} | awk 'NR % 4 == 2' > {output}",
    )
    args = parser.parse_args()
    return args


def sam_fastq_tee(input, tag_out, header_out, fq_out, fmt="fastq"):
    from mrfifo.plumbing import open_named_pipe

    n = 0
    # print("opening input", input, type(input))
    input = open_named_pipe(input, "rt")
    # print("opening header_out")
    header_out = open_named_pipe(header_out, "wt")
    # print("opening tag_out")
    tag_out = open_named_pipe(tag_out, "wt")
    # print("opening fq_out")
    fq_out = open_named_pipe(fq_out, "wt")
    # print("done. reading sam...")
    for sam in input:
        # print(sam)
        if sam.startswith("@"):
            header_out.write(sam)
            continue

        elif header_out is not None:
            header_out.flush()
            header_out.close()
            header_out = None

        parts = sam.rstrip().split("\t")
        qname = parts[0]
        seq = parts[9]
        qual = parts[10]
        tags = parts[11:]

        tag_out.write("\t".join(tags) + "\n")
        # fq_out.write(f"@{qname}\n{seq}\n+\n{qual}\n")
        if fmt == "fastq":
            fq_out.write(f"@{qname}\n{seq}\n+\n{qual}\n")
        else:
            fq_out.write(f">{qname}\n{seq}\n")
        n += 1

    fq_out.flush()
    fq_out.close()

    tag_out.flush()
    tag_out.close()

    return n


def wrap_cmd(input, output, cmd="cat {input} | tr '@' '\t' > {output}"):
    import os

    os.system(cmd.format(input=input, output=output))


def merge_output(original, aligned, header):
    # print("reading header")
    for line in header:
        print(line.rstrip())

    # print("done")

    for ori, ali in zip(original, aligned):
        print(ali.rstrip() + "\t" + ori.rstrip())


def print_stream(input, outpath):
    from mrfifo.plumbing import open_named_pipe

    input = open_named_pipe(input, "rt")
    out = open(outpath, "w")

    for line in input:
        out.write(line)

    input.flush()
    input.close()

    out.flush()
    out.close()


# the order in which the fifos are opened for reading or writing matters if more than r or w each are opened
# they need to be opened for reading in the same order that they are opened for writing. Else you can get a deadlock upon initialization.


def main(args):
    w = (
        mf.Workflow("bamwrap", total_pipe_buffer_MB=8)
        .BAM_reader(input=args.bam_in)
        .add_job(
            func=sam_fastq_tee,
            input=mf.FIFO("input_sam", "r"),
            tag_out=mf.FIFO("sam_tags", "w"),
            header_out=mf.FIFO("header_out", "w"),
            fq_out=mf.FIFO("fastq_input", "w"),
            fmt=args.format,
            _manage_fifos=False,
            job_name="{workflow}.tee{n}",
        )
        .add_job(
            func=wrap_cmd,
            input=mf.FIFO("fastq_input", "r"),
            output=mf.FIFO("sam_output", "w"),
            cmd=args.cmd,
            _manage_fifos=False,
            job_name="{workflow}.wrap_cmd{n}",
        )
        .add_job(
            func=merge_output,
            header=mf.FIFO("header_out", "r"),
            original=mf.FIFO("sam_tags", "r"),
            aligned=mf.FIFO("sam_output", "r"),
        )
    ).run()

    return w.result_dict


def cmdline():
    args = parse_args()
    import spacemake.util as util

    util.setup_logging(args, name="spacemake.bin.bamwrap")
    return main(args)


if __name__ == "__main__":
    cmdline()
