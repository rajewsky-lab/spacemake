from spacemake.contrib import __version__, __license__, __author__, __email__
import spacemake.util as util

import logging


def load_fasta(fname):
    logging.info(f"loading fasta '{fname}'")
    chrom_sizes = {}
    seq_d = {}
    full_names = {}
    for name, seq in util.fasta_chunks(open(fname)):
        fa_id = name.split()[0]
        full_names[fa_id] = name
        seq_d[fa_id] = seq.upper()
        chrom_sizes[fa_id] = len(seq)

    return seq_d, chrom_sizes, full_names


def match_exp(regexp, s, no_match="na"):
    if not "(" in regexp:
        return regexp

    import re
    M = re.search(regexp, s)
    if not M:
        return no_match
    else:
        return M.groups()[0]


def render_gff(
    chrom,
    start,
    end,
    strand="+",
    score=".",
    type="exon",
    source="make_gtf_from_fa.py",
    phase=".",
    **kw,
):
    attrs = [f'{key} "{value}";' for key, value in sorted(kw.items())]
    attr_str = " ".join(attrs)
    columns = [
        chrom,
        source,
        type,
        str(start + 1),
        str(end),
        str(score),
        strand,
        str(phase),
        attr_str,
    ]
    return "\t".join(columns)


def parse_args():
    parser = util.make_minimal_parser()
    parser.add_argument("fasta", help="fasta-input file")
    parser.add_argument("--source", help="value for source column", default="make_gtf_from_fa.py")
    parser.add_argument("--gene_name", default=r"(^\S+)", help="regexp to match against fasta_id to extract gene_name value")
    parser.add_argument("--transcript_id", default=r"(^\S+)", help="regexp to match against fasta_id to extract transcript_id")
    parser.add_argument("--transcript_type", default=r"na", help="regexp to match against fasta_id to extract transcript_type value (default='na')")
    parser.add_argument("--gtf", default="/dev/stdout", help="output GTF here (default=stdout)")
    parser.add_argument("--features", default="exon,transcript", help="comma-separated list of features to generate. default='exon,transcript'")
    args = parser.parse_args()
    return args


def main(args):
    seq_d, chrom_sizes, full_names = load_fasta(args.fasta)
    with open(args.gtf, "w") as out:
        for chrom, size in sorted(chrom_sizes.items()):
            name = full_names[chrom]
            for t in args.features.split(','):
                gff = render_gff(
                    chrom,
                    0,
                    size,
                    source=args.source,
                    type=t,
                    gene_id=chrom,
                    gene_name=match_exp(args.gene_name, name),
                    transcript_id=match_exp(args.transcript_id, name),
                    transcript_type=match_exp(args.transcript_type, name, no_match="other"),
                    # transcript_name=chrom,
                )
                out.write(gff + '\n')

if __name__ == "__main__":
    args = parse_args()
    util.setup_logging(args, name="spacemake.bin.make_gtf_from_fa")
    main(args)
    