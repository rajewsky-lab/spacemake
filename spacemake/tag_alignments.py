import numpy as np
import logging
import pandas as pd

# import cutadapt.align
from collections import defaultdict

__version__ = "0.9"
__author__ = ["Marvin Jens"]
__license__ = "GPL"


def parse_cmdline():
    import argparse

    parser = argparse.ArgumentParser(
        description="tag alignments in a BAM file using various rules (aRNA and miRNA index)"
    )
    parser.add_argument(
        "bam_in",
        help="bam input (default=stdin)",
        default="/dev/stdin",
        # nargs="+",
    )
    parser.add_argument(
        "--bam-out",
        help="bam output (default=stdout)",
        default="/dev/stdout",
    )
    parser.add_argument(
        "--bam-out-mode",
        help="bam output mode (default=b)",
        default="b",
    )
    parser.add_argument(
        "--rules",
        help="which rules to apply. Options are 'aRNA','miRNA'",
        default="aRNA",
        choices=["aRNA", "miRNA", "genome"],
    )
    parser.add_argument(
        "--lookup",
        help="tab-separated table with reference -> function, name mappings",
        default="",
    )
    parser.add_argument(
        "--primer-pool",
        help="tab-separated table with reference -> function, name mappings",
        default="",
    )

    parser.add_argument(
        "--skim",
        help="skim through the BAM by investigating only every <skim>-th record (default=1 off)",
        default=1,
        type=int,
    )
    # parser.add_argument(
    #     "--min-length",
    #     help="minimal allowed read-length left after trimming (default=18)",
    #     type=int,
    #     default=18,
    # )
    # parser.add_argument(
    #     "--min-qual",
    #     help="minimal quality score for quality trimming (default=30)",
    #     type=int,
    #     default=20,
    # )
    parser.add_argument(
        "--stats-out",
        help="write tab-separated table with tagging overview results here",
        default="",
    )
    return parser.parse_args()


def make_header(bam):
    import os
    import sys

    header = bam.header.to_dict()
    progname = os.path.basename(__file__)
    # if "PG" in header:
    # for pg in header['PG']:
    #     if pg["ID"] == progname:
    #         progname = progname + ".1"

    pg_list = header.get("PG", [])
    pg = {
        "ID": progname,
        "PN": progname,
        "CL": " ".join(sys.argv[1:]),
        "VN": __version__,
    }
    if len(pg_list):
        pg["PP"] = pg_list[-1]["ID"]

    header["PG"] = pg_list + [pg]
    return header


class AbundantRNATagger:
    def __init__(self, bam, lkup_table={}):
        self.logger = logging.getLogger("AbundantRNATagger")
        self.lkup_table = lkup_table
        self.bam = bam
        self.tid_lkup = {}
        self.umi = set()

        self.counter = defaultdict(int)
        self.tag_names_to_count = ["af", "an"]

    def cached_name_for_tid(self, tid):
        if not tid in self.tid_lkup:
            self.tid_lkup[tid] = self.bam.get_reference_name(tid)
        return self.tid_lkup[tid]

    def make_tags(self, aln):
        name = self.cached_name_for_tid(aln.tid)
        _tags = aln.get_tags()
        tags = dict(_tags)
        key = (aln.query_sequence, tags.get("CB", "NA"), tags.get("MI", "NA"))

        if key in self.umi:
            ax = "PCR"
        else:
            ax = "UMI"

        self.umi.add(key)
        af, an = self.lkup_table.get(name, ("NA", "NA"))
        new_tags = [("ax", ax), ("an", an), ("af", af)]
        tags.update(dict(new_tags))

        return name, _tags + new_tags, tags

    def count(self, tags):
        self.counter[("reads", "stats", "N_reads")] += 1
        for n in self.tag_names_to_count:
            self.counter[("reads", n, tags[n])] += 1

        if tags["ax"] == "UMI":
            self.counter[("UMIs", "stats", "N_UMIs")] += 1
            for n in self.tag_names_to_count:
                self.counter[("UMIs", n, tags[n])] += 1

    def tag_alignment(self, aln):
        name, _tags, tags = self.make_tags(aln)
        self.count(tags)
        aln.set_tags(_tags)
        return aln, tags

    def write_stats(self, fname):
        N_reads = self.counter[("reads", "stats", "N_reads")]
        N_UMIs = self.counter[("UMIs", "stats", "N_UMIs")]

        with open(fname, "wt") as f:
            f.write("scope\ttag\tname\tcount\tpercent\n")
            for (scope, tag, name), v in sorted(
                self.counter.items(), key=lambda x: (x[0][:2], -x[1])
            ):
                if scope == "reads":
                    total = N_reads
                else:
                    total = N_UMIs
                f.write(f"{scope}\t{tag}\t{name}\t{v}\t{100.0 * v/total:.2f}\n")


class miRNATagger(AbundantRNATagger):
    def __init__(self, bam, lkup_table, lkup_targeted):
        AbundantRNATagger.__init__(self, bam, lkup_table)

        self.lkup_targeted = lkup_targeted
        self.lkup_targeted_gene = {}
        for key, pools in lkup_targeted.items():
            gene = self.lkup_table[key]
            self.lkup_targeted_gene[gene] = pools

        self.tag_names_to_count.append("aa")

    def make_tags(self, aln):
        name, _tags, tags = AbundantRNATagger.make_tags(self, aln)

        an = tags["an"]
        # parse CIGAR to get 5'/3' clipping and longest contiguous match
        n_match = 0
        n_clip5 = 0
        n_clip3 = 0
        if aln.cigartuples:
            last_i = len(aln.cigartuples) - 1
            for i, (op, n) in enumerate(aln.cigartuples):
                if op == 0:
                    n_match = max(n, n_match)

                if (i == 0) and (op == 4):
                    n_clip5 = n

                if (i == last_i) and (op == 4):
                    n_clip3 = n

        # populate "artifact warning" tag
        aa = []
        if "TSO" in tags.get("A5", ""):
            aa.append("TSO")
            if name in self.lkup_targeted:
                # aa.append(f"POOL:{self.lkup_targeted[name]}")
                aa.append(f"POOL")

            elif an in self.lkup_targeted_gene:
                aa.append(f"pool")

        if tags["af"] == "junk":
            aa.append("junk")

        if not "polyA" in tags.get("A3", ""):
            aa.append("no_polyA")

        if n_clip5 != 2:
            aa.append("no_lig_NN")

        if n_match < 17:
            aa.append("short_match")

        if n_clip3 > 10:
            aa.append("too_much_clipped")

        if not aa:
            aa.append("pass")

        aa = ",".join(aa)

        new_tags = [
            ("aa", aa),
            ("c5", n_clip5),
            ("c3", n_clip3),
            ("cm", n_match),
        ]

        tags.update(dict(new_tags))
        return (name, _tags + new_tags, tags)


class GenomeTagger(AbundantRNATagger):
    def select_genes(self, strand, gn_list, gf_list, gs_list):
        prio = {
            # sense-match, exonic
            (True, True): 0,
            (True, False): 1,
            (False, True): 20,
            (False, False): 30,
        }
        results = []
        best = 10000
        for i, (gn, gf, gs) in enumerate(zip(gn_list, gf_list, gs_list)):
            sense_match = gs == strand
            exonic = (gf == "UTR") or (gf == "CODING")
            p = prio[(sense_match, exonic)]
            results.append((p, gn, i))
            best = min(best, p)

        genes = set()
        gf = ""
        for p, gn, i in results:
            if p == best:
                genes.add(gn)
                gf = gf_list[i]

        genes = list(genes)
        if len(genes) == 1:
            return genes[0], gf
        else:
            return ",".join(sorted(genes)), gf

    def make_tags(self, aln):
        # name = aln.get_tag("gn")
        _tags = aln.get_tags()
        tags = dict(_tags)

        name, af = self.select_genes(
            "-" if aln.is_reverse else "+",
            tags.get("gn", "").split(","),
            tags.get("gf", "").split(","),
            tags.get("gs", "").split(","),
        )

        key = (aln.query_sequence[:20], tags.get("CB", "NA"), tags.get("MI", "NA"))

        if key in self.umi:
            ax = "PCR"
        else:
            ax = "UMI"

        self.umi.add(key)

        new_tags = [("ax", ax), ("an", name), ("af", af)]
        tags.update(dict(new_tags))

        return name, _tags + new_tags, tags


if __name__ == "__main__":
    args = parse_cmdline()

    stats = defaultdict(int)
    total = defaultdict(int)
    lhist = defaultdict(int)

    import pysam

    bam_in = pysam.AlignmentFile(args.bam_in, "rb", check_sq=False)
    bam_out = pysam.AlignmentFile(
        args.bam_out, f"w{args.bam_out_mode}", header=make_header(bam_in)
    )

    lkup_table = {}
    if args.lookup:
        for row in pd.read_csv(args.lookup, sep="\t").itertuples():
            lkup_table[row.ref] = (row.af, row.an)

    lkup_pool = {}
    if args.primer_pool:
        for row in pd.read_csv(args.primer_pool, sep="\t").itertuples():
            txt = []
            # print(row)
            if row._2:
                txt.append("3")
            if row._3:
                txt.append("12")
            if row._4:
                txt.append("22")
            if row._5:
                txt.append("mmu")

            lkup_pool[row.name] = ",".join(txt)

    if args.rules == "aRNA":
        tagger = AbundantRNATagger(bam_in, lkup_table)
    elif args.rules == "miRNA":
        tagger = miRNATagger(bam_in, lkup_table, lkup_pool)
    else:
        tagger = GenomeTagger(bam_in)

    for i, aln in enumerate(bam_in.fetch(until_eof=True)):
        if args.skim and i % args.skim != 0:
            continue

        aln, tags = tagger.tag_alignment(aln)
        bam_out.write(aln)

    if args.stats_out:
        tagger.write_stats(args.stats_out)
