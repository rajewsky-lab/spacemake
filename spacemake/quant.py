import pysam
import numpy as np
import sys
import logging
import argparse
from collections import defaultdict


def out_counts_bulk(f, counts, discard, stats):
    discard_keys = sorted(discard.keys())
    f.write("name\tL\tcount\tdropped\t" + "\t".join(discard_keys) + "\n")
    all_keys = stats.keys()
    all_counts = dict([(k, counts[k]) for k in all_keys])
    for mir, count in sorted(all_counts.items(), key=lambda x: -x[1]):
        # L =
        # print(f"{mir}\t{sam.get_reference_length(mir)}\t{count}\t{rev[mir]}")
        f.write(
            f"{mir}\t-1\t{count}\t{stats[mir]['n_dropped']}\t"
            + "\t".join([str(discard[reason][mir]) for reason in discard_keys])
            + "\n"
        )


class AlignmentClassifier:
    def __init__(
        self,
        sample_name,
        bam_in,
        gene_assign_mode="chrom",
        chrom_to_gene={},
        ignore_gf=["INTRONIC", "INTERGENIC"],
    ):
        self.logger = logging.getLogger(f"AlignmentClassifier({sample_name})")
        self.sample_name = sample_name
        self.bam = pysam.AlignmentFile(bam_in, check_sq=False)
        self.gene_names = {}
        self.chrom_to_gene = chrom_to_gene
        self.ignore_gf = set(ignore_gf)
        self.get_gene = {
            "chrom": self.fast_refname,
            "gn_tag": self.parse_gn_gf,
        }[gene_assign_mode]

        self.uniq_reads = set()

    def fast_refname(self, aln):
        if not aln.tid in self.gene_names:
            chrom = self.bam.get_reference_name(aln.tid)
            gene = self.chrom_to_gene.get(chrom, chrom)
            self.gene_names[aln.tid] = gene

        return self.gene_names[aln.tid]

    def parse_gn_gf(self, aln):
        if not aln.has_tag("gf"):
            return "NA"

        genes = set()
        for gf, gn in zip(aln.get_tag("gf").split(","), aln.get_tag("gn").split(",")):
            if gf not in self.ignore_gf:
                genes.add(gn)

        genes = list(genes)
        if len(genes) == 1:
            return genes[0]
        else:
            return "multi_mapper"

    def is_uniq(self, aln):
        key = (aln.query_sequence, aln.get_tag("CB"), aln.get_tag("MI"))
        u = key in self.uniq_reads
        self.uniq_reads.add(u)
        return u

    def __iter__(self):
        for aln in self.bam.fetch(until_eof=True):
            yield ClassifiedAlignment(
                self, aln, self.get_gene(aln), is_dup=not self.is_uniq(aln)
            )


class ClassifiedAlignment:
    def __init__(self, parent, aln, gene, is_dup=False):
        self.aln = aln
        self.parent = parent
        self.tags = dict(aln.get_tags())
        self.flags = set()
        self.cell = self.tags.get("CB", "NA")
        self.umi = self.tags.get("MI", "NA")
        self.gene = gene
        self.is_dup = is_dup

    def check_qname(self, keywords=[]):
        for n in keywords:
            if n in self.aln.qname:
                self.flags.add(f"name_has_{n}")
            else:
                self.flags.add(f"name_missing_{n}")

        return self

    def check_tags(self, **kw):
        for tag, keywords in kw.items():
            tagval = self.tags.get(tag, "")
            for n in keywords:
                if n in tagval:
                    self.flags.add(f"{tag}_has_{n}")
                else:
                    self.flags.add(f"{tag}_missing_{n}")

        return self

    def check_CIGAR(self):
        n_match = 0
        c5 = 0
        c3 = 0
        cigar = self.aln.cigartuples
        if cigar:
            op, n = cigar[0]
            if op == 4:  # soft-clip
                c5 = n
            elif op == 0:  # match
                n_match = max(n_match, n)

        if len(cigar) > 1:
            for op, n in cigar[1:]:
                if op == 4:  # soft-clip again, must be 3' end
                    c3 = n
                elif op == 0:  # match
                    n_match = max(n_match, n)

        self.clip5 = c5
        self.clip3 = c3
        self.n_match = n_match

        return self


class DGE:
    def __init__(self, _assert=False, main_channel="count"):
        self.main_channel = main_channel
        self.DGE_umis = defaultdict(lambda: defaultdict(set))
        self.DGE_reads = defaultdict(lambda: defaultdict(int))
        self.DGE_cells = set()
        self.DGE_genes = set()

        # marginal counters for each cell
        self.DGE_cell_reads = defaultdict(int)
        self.DGE_cell_umis = defaultdict(set)
        self.DGE_cell_genes = defaultdict(set)

        # marginal counts for each gene
        self.DGE_gene_reads = defaultdict(int)
        self.DGE_gene_umis = defaultdict(set)
        self.DGE_gene_cells = defaultdict(set)

        # flag to enable consistency checking (for debugging)
        self._assert = _assert
        self.channels = set([main_channel])

    def add_read(self, gene, cell, umi, channel="count"):
        self.DGE_cells.add(cell)
        self.DGE_genes.add(gene)

        # margin counters for validation
        if channel == self.main_channel:
            self.DGE_cell_reads[cell] += 1

        if (gene, umi) in self.DGE_cell_umis[cell]:
            dup = True
        else:
            dup = False

        self.DGE_cell_umis[cell].add((gene, umi))
        if self._assert:
            self.DGE_cell_genes[cell].add(gene)
            self.DGE_gene_reads[gene] += 1
            self.DGE_gene_umis[gene].add((cell, umi))
            self.DGE_gene_cells[gene].add(cell)

        # main count matrices
        self.DGE_umis[(gene, cell)][channel].add(umi)
        # if channel == self.main_channel:
        self.DGE_reads[(gene, cell)][channel] += 1

        self.channels.add(channel)

        return dup

    def make_DGEs(self):
        import scipy.sparse
        import anndata

        # count UMIs, reads across all channels into sparse arrays
        ind_of_cell = {}
        obs = []
        for cell in sorted(self.DGE_cells):
            ind_of_cell[cell] = len(obs)
            obs.append(cell)

        ind_of_gene = {}
        var = []
        for gene in sorted(self.DGE_genes):
            ind_of_gene[gene] = len(var)
            var.append(gene)

        umi_counts_by_channel = defaultdict(list)
        read_counts_by_channel = defaultdict(list)
        row_ind = []
        col_ind = []
        for (gene, cell), umi_dict in self.DGE_umis.items():
            for channel in self.channels:
                umis = umi_dict[channel]
                umi_counts_by_channel[channel].append(
                    len(umis)
                )  # <- the actual expression count
                read_counts_by_channel[channel].append(
                    self.DGE_reads[(gene, cell)][channel]
                )

            row_ind.append(ind_of_cell[cell])
            col_ind.append(ind_of_gene[gene])

        # print(
        #     ">>>DIMENSION CHECK",
        #     self.channels,
        #     len(self.DGE_umis.items()),
        #     [(k, len(v)) for k, v in umi_counts_by_channel.items()],
        #     len(umi_counts_by_channel[self.main_channel]),
        #     len(row_ind),
        #     len(col_ind),
        # )
        X = scipy.sparse.csr_matrix(
            (umi_counts_by_channel[self.main_channel], (row_ind, col_ind))
        )
        # convert to anndata object and store
        adata = anndata.AnnData(X)
        adata.obs_names = obs
        adata.var_names = var
        adata.layers[f"reads_{self.main_channel}"] = scipy.sparse.csr_matrix(
            (read_counts_by_channel[self.main_channel], (row_ind, col_ind))
        )

        for channel in self.channels:
            if channel == self.main_channel:
                continue

            umi_counts = umi_counts_by_channel[channel]

            adata.layers[channel] = scipy.sparse.csr_matrix(
                (umi_counts, (row_ind, col_ind))
            )
            read_counts = read_counts_by_channel[channel]
            adata.layers[f"reads_{channel}"] = scipy.sparse.csr_matrix(
                (read_counts, (row_ind, col_ind))
            )
            # Do we want to keep the extra UMI or read counts?
            # adata.layers[name] = aextra.layers["reads"]

        # print(adata)
        return adata

    def check_DGEs_vs_margin_counts(self, ann_umis, ann_reads):
        import numpy as np

        if not self._assert:
            raise ValueError(
                "DGE was not initialized with _assert = True, which is required for check_DGEs_vs_margin_counts()"
            )

        cell_names = ann_umis.obs_names
        gene_names = ann_umis.var_names

        def sets_to_vec(S, names):
            return np.array([len(S[n]) for n in names])

        def count_to_vec(S, names):
            return np.array([S[n] for n in names])

        def flatvec(m):
            return np.array(m).ravel()

        assert (
            flatvec(ann_umis.X.sum(axis=1))
            == sets_to_vec(self.DGE_cell_umis, cell_names)
        ).all()
        assert (
            flatvec(ann_umis.X.sum(axis=0))
            == sets_to_vec(self.DGE_gene_umis, gene_names)
        ).all()

        assert (
            flatvec(ann_reads.X.sum(axis=1))
            == count_to_vec(self.DGE_cell_reads, cell_names)
        ).all()
        assert (
            flatvec(ann_reads.X.sum(axis=0))
            == count_to_vec(self.DGE_gene_reads, gene_names)
        ).all()

        detect_reads = ann_reads.X > 0
        assert (
            flatvec(detect_reads.sum(axis=1))
            == sets_to_vec(self.DGE_cell_genes, cell_names)
        ).all()
        assert (
            flatvec(detect_reads.sum(axis=0))
            == sets_to_vec(self.DGE_gene_cells, gene_names)
        ).all()

        detect_umis = ann_umis.X > 0
        assert (
            flatvec(detect_umis.sum(axis=1))
            == sets_to_vec(self.DGE_cell_genes, cell_names)
        ).all()
        assert (
            flatvec(detect_umis.sum(axis=0))
            == sets_to_vec(self.DGE_gene_cells, gene_names)
        ).all()


# def assess_reads_per_UMI(ann_umis, ann_reads):
#     rpu = (
#         np.array(ann_reads.X.sum(axis=1), dtype=float).ravel()
#         / np.array(ann_umis.X.sum(axis=1), dtype=float).ravel()
#     )
#     perc = np.percentile(rpu, [5, 25, 50, 75, 95])
#     # print(perc)
#     return perc


def parse_cmdline():
    parser = argparse.ArgumentParser(
        description="quantify miRNA expression from BAM files"
    )

    parser.add_argument(
        "bam_in",
        help="bam input (default=stdin)",
        default="/dev/stdin",
        # nargs="+",
    )
    parser.add_argument(
        "--sample-name",
        help="sample identifier (default='NA')",
        default="NA",
    )

    parser.add_argument(
        "--parse-gn",
        help="parse gf,gn BAM tags instead of mapping information",
        default=False,
        action="store_true"
        # nargs="+",
    )
    parser.add_argument(
        "--ignore-gf",
        help="which gene functions to ignore? comma-separated string. default=INTRONIC,INTERGENIC",
        default="INTRONIC,INTERGENIC",
        # nargs="+",
    )

    parser.add_argument(
        "--translate",
        help="translate gene-names using this two-column lookup table w columns 'original' and 'target'",
        default="",
        # nargs="+",
    )
    parser.add_argument(
        "--skim",
        help="skim through the BAM by investigating only every <skim>-th record (default=1 off)",
        default=1,
        type=int
        # nargs="+",
    )
    parser.add_argument(
        "--output-counted",
        help="direct neg control alignments that are counted to this BAM file. default=off",
        default="",
    )
    parser.add_argument(
        "--output-discarded",
        help="direct reads that are discarded for reasons other than PCR-duplicate to this BAM file. default=off",
        default="",
    )
    parser.add_argument(
        "--output-passed",
        help="direct reads that pass all filters (but may be PCR-duplicates) to this BAM file. default=off",
        default="",
    )

    parser.add_argument(
        "--output-detailed",
        help="direct counts for individual miRNA family members to this file default=off",
        default="",
    )
    parser.add_argument(
        "--output-DGE",
        help="output a single-cell digital gene expression matrix with UMI COUNTS",
        default="",
    )
    parser.add_argument(
        "--layers",
        default="reads",
        help="comma-separated list of extra info to add as layers",
    )
    # parser.add_argument(
    #     "--output-DGE-reads",
    #     help="output a single-cell digital gene expression matrix with READ COUNTS",
    #     default="",
    # )
    parser.add_argument(
        "--count-func",
        default="ligation_product",
        help="the function which decides whether to count a read",
        choices=["ligation_product", "targeted_primer", "everything"],
    )
    parser.add_argument(
        "--min-match",
        default=17,
        type=int,
        help="minimum number of matching bases required (default=17)",
    )
    parser.add_argument(
        "--max-clipped",
        default=50,
        type=int,
        help="maximum number of soft-clipped 3'bases allowed (default=50)",
    )
    parser.add_argument(
        "--min-BC",
        default=0,
        type=int,
        help="expect at least this number of barcode bases at the 5' end (default=0)",
    )
    parser.add_argument(
        "--name-filter",
        default="",
        help="require presence of this string in qname (default='' -> off)",
    )
    parser.add_argument(
        "--name-exclude",
        default="",
        help="discard if this string is present in qname (default='' -> off)",
    )

    return parser.parse_args()


def count_ligation_product(ca):
    return (
        not (ca.short or ca.primer) and "A3_has_polyA" in ca.flags and ca.gene != "NA"
    )


def count_targeted_primer(ca):
    return ca.primer


def count_everything(ca):
    return True


if __name__ == "__main__":
    args = parse_cmdline()

    lkup = {}
    if args.translate:
        for row in pd.read_csv(args.translate, sep="\t").itertuples():
            lkup[row.original] = row.target

    if args.parse_gn:
        gene_mode = "gn_tag"
    else:
        gene_mode = "chrom"

    ignore = args.ignore_gf.split(",")

    count_func = {
        "ligation_product": count_ligation_product,
        "targeted_primer": count_targeted_primer,
        "everything": count_everything,
    }[args.count_func]
    dge = DGE()
    for ca in AlignmentClassifier(
        args.sample_name,
        args.bam_in,
        chrom_to_gene=lkup,
        gene_assign_mode=gene_mode,
        ignore_gf=ignore,
    ):
        ca = (
            # .check_qname(keywords=["TSO", "polyA"])
            ca.check_tags(
                A3=["polyA"], A5=["TSO"]
            ).check_CIGAR()  # fills flags["A3_has_polyA"] etc.  # fills .clip5, .n_match, .clip3
        )

        ca.short = ca.n_match < args.min_match
        ca.reverse = ca.aln.is_reverse
        # ca.PCR_dup  = ca.
        ca.primer = ("A5_has_TSO" in ca.flags) or ("A3_missing_polyA" in ca.flags)
        ca.count = count_func(ca)
        # print(ca.n_match, "flags=", ca.flags)
        # print(ca.short, ca.primer, ca.reverse, ca.count)

        channels = ["count", "short", "reverse", "primer"]
        for c in channels:
            if getattr(ca, c):
                # print(ca.aln)
                # print("counting as", c)
                dge.add_read(gene=ca.gene, cell=ca.cell, umi=ca.umi, channel=c)

        # print("next")
    if args.output_DGE and len(dge.DGE_reads):
        adata = dge.make_DGEs()
        # dge.check_DGEs_vs_margin_counts(ann_umis, ann_reads)
        # print("storing AnnData object")
        adata.write(args.output_DGE)
        # ann_reads.write(args.output_DGE_reads)

        # sys.stderr.write(
        #     f"### reads-to-UMI ratio quartiles: {assess_reads_per_UMI(ann_umis, ann_reads)} \n"
        # )

    # # print("starting")
    # counts = defaultdict(int)
    # discarded = defaultdict(lambda: defaultdict(int))  # key is [reason][miRNA] -> count
    # counts_detailed = defaultdict(int)
    # discarded_detailed = defaultdict(lambda: defaultdict(int))
    # rev = defaultdict(int)
    # rev_detailed = defaultdict(int)
    # stats = defaultdict(lambda: defaultdict(int))
    # stats_detailed = defaultdict(lambda: defaultdict(int))
    # clipped = defaultdict(lambda: defaultdict(int))
    # nclipped = defaultdict(int)
    # clipped5 = defaultdict(lambda: defaultdict(int))
    # nclipped5 = defaultdict(int)

    # sam = pysam.AlignmentFile(args.bam_in, "rb", check_sq=False)
    # if args.output_counted:
    #     sam_counted = pysam.AlignmentFile(args.output_counted, "wb", template=sam)
    # else:
    #     sam_counted = None

    # if args.output_discarded:
    #     sam_discard = pysam.AlignmentFile(args.output_discarded, "wb", template=sam)
    # else:
    #     sam_discard = None

    # if args.output_passed:
    #     sam_passed = pysam.AlignmentFile(args.output_passed, "wb", template=sam)
    # else:
    #     sam_passed = None

    # full_names = set()
    # simple_names = set()
    # dge = DGE()
    # layer_names = set(args.layers.split(","))
    # # print(f"DGE-counting the following events: {layer_names}")

    # def record_dge_extra(full, reason, read):
    #     if reason in layer_names:
    #         dge.add_read(full, read.get_tag("CB"), read.get_tag("MI"), channel=reason)

    # def flag_drop(simple, full, reason, read):
    #     stats[simple][f"n_{reason}"] += 1
    #     stats_detailed[full][f"n_{reason}"] += 1
    #     discarded[reason][simple] += 1
    #     discarded_detailed[reason][full] += 1

    #     record_dge_extra(full, reason, read)
    #     record_dge_extra(full, "dropped", read)

    #     if reason != "PCR_duplicate" and sam_discard:
    #         key = read.query_sequence + read.get_tag("CB") + read.get_tag("MI")
    #         # if not key in discarded_uniq:
    #         read.set_tag("XD", reason)
    #         sam_discard.write(read)

    #         # discarded_uniq.add(key)

    #     return True

    # for i_read, read in enumerate(sam.fetch(until_eof=True)):
    #     if (i_read % args.skim) != 0:
    #         continue

    #     # print(read)
    #     full, simple, L = get_mirname(sam, read, args.parse_gf)
    #     full_names.add(full)
    #     simple_names.add(simple)
    #     stats[simple]["n_reads"] += 1
    #     stats_detailed[full]["n_reads"] += 1
    #     if full != simple:
    #         stats[full]["n_reads"] += 1
    #     drop = False

    #     n_match = 0
    #     if args.min_match:
    #         for op, n in read.cigartuples:
    #             # find the largest contiguous stretch of aligned bases
    #             if op == 0:
    #                 n_match = max(n_match, n)

    #     if n_match < args.min_match:
    #         drop = flag_drop(simple, full, "match_too_short", read)

    #     n5 = 0
    #     if read.cigartuples:
    #         op, n = read.cigartuples[0]
    #         if op == 4:
    #             n5 = n

    #     if "TSO" in read.qname:
    #         if full != simple:
    #             clipped5[full]["TATGGG"] += 1
    #         clipped5[simple]["TATGGG"] += 1
    #         if full != simple:
    #             nclipped5[full] += 1
    #         nclipped5[simple] += 1
    #     elif n5:
    #         n_clipped = n5
    #         if full != simple:
    #             clipped5[full][read.query_sequence[:n]] += 1
    #         clipped5[simple][read.query_sequence[:n]] += 1
    #         if full != simple:
    #             nclipped5[full] += 1
    #         nclipped5[simple] += 1

    #     n_clipped = 0
    #     if read.cigartuples:
    #         op, n = read.cigartuples[-1]
    #         if op == 4:
    #             clipped_seq = clip_pA(read.query_sequence[-n:])
    #             n_clipped = len(clipped_seq)
    #             if full != simple:
    #                 clipped[full][clipped_seq] += 1
    #             clipped[simple][clipped_seq] += 1
    #         else:
    #             if full != simple:
    #                 clipped[full][""] += 1
    #             clipped[simple][""] += 1

    #     if full != simple:
    #         nclipped[full] += 1
    #     nclipped[simple] += 1
    #     clipped[simple]["all_occurrences"] += 1
    #     if full != simple:
    #         clipped[full]["all_occurrences"] += 1

    #     if n_clipped > args.max_clipped:
    #         drop = flag_drop(simple, full, "too_much_clipped", read)

    #     if args.min_BC and (n5 != args.min_BC):
    #         drop = flag_drop(simple, full, "missing_BC", read)

    #     if args.name_filter:
    #         if not args.name_filter in read.qname:
    #             drop = flag_drop(simple, full, f"name_missing_{args.name_filter}", read)
    #             # stats[simple]["n_name_filtered"] += 1

    #             # stats_detailed[full]["n_name_filtered"] += 1
    #             # discarded[f"name_missing_{args.name_filter}"][simple] += 1
    #             # discarded_detailed[f"name_missing_{args.name_filter}"][full] += 1

    #             # drop = True

    #     if args.name_exclude:
    #         if args.name_exclude in read.qname:
    #             drop = flag_drop(simple, full, f"name_has_{args.name_exclude}", read)

    #             # stats[simple]["n_name_excluded"] += 1
    #             # stats_detailed[full]["n_name_excluded"] += 1
    #             # discarded[f"name_has_{args.name_exclude}"][simple] += 1
    #             # discarded_detailed[f"name_has_{args.name_exclude}"][full] += 1

    #             # drop = True

    #     if read.is_reverse:
    #         drop = flag_drop(simple, full, "reverse", read)

    #         # stats[simple]["n_reverse"] += 1
    #         # stats_detailed[full]["n_reverse"] += 1
    #         # discarded["reverse"][simple] += 1
    #         # discarded_detailed["reverse"][full] += 1

    #     if not drop:
    #         # dup = dge.add_read(full, read.get_tag("CB"), read.get_tag("MI"))
    #         # if dup:
    #         #     drop = flag_drop(simple, full, "PCR_duplicate")
    #         if sam_passed:
    #             sam_passed.write(read)

    #         if args.output_DGE:
    #             dup = dge.add_read(full, read.get_tag("CB"), read.get_tag("MI"))
    #             if dup:
    #                 drop = flag_drop(simple, full, "PCR_duplicate", read)

    #     if drop:
    #         stats[simple]["n_dropped"] += 1
    #         stats_detailed[full]["n_dropped"] += 1
    #     else:
    #         counts[simple] += 1
    #         counts_detailed[full] += 1
    #         if simple.startswith("NegCtrl") and sam_counted:
    #             sam_counted.write(read)

    #         stats[simple]["n_kept"] += 1
    #         stats_detailed[full]["n_kept"] += 1

    # if args.output_DGE and len(dge.DGE_reads):
    #     adata = dge.make_DGEs()
    #     # dge.check_DGEs_vs_margin_counts(ann_umis, ann_reads)
    #     # print("storing AnnData object")
    #     adata.write(args.output_DGE)
    #     # ann_reads.write(args.output_DGE_reads)

    #     # sys.stderr.write(
    #     #     f"### reads-to-UMI ratio quartiles: {assess_reads_per_UMI(ann_umis, ann_reads)} \n"
    #     # )

    # out_counts_bulk(sys.stdout, counts, discarded, stats)
    # # for mir, count in sorted(counts.items(), key=lambda x: -x[1]):
    # #     # print(f"{mir}\t{sam.get_reference_length(mir)}\t{count}\t{rev[mir]}")
    # #     print(
    # #         f"{mir}\t-1\t{count}\t{rev[mir]}\t{stats[mir]['n_dropped']}\t{stats[mir]['n_dropped']/count:.1f}"
    # #     )

    # if args.output_detailed:
    #     out_counts_bulk(
    #         open(args.output_detailed, "wt"),
    #         counts_detailed,
    #         discarded_detailed,
    #         stats_detailed,
    #     )

    # for simple_name in simple_names:
    #     data = stats[simple_name]
    #     for key, value in sorted(data.items()):
    #         perc = 100 * value / data["n_reads"]
    #         sys.stderr.write(f"{simple_name}\t{key}\t{value}\t{perc:.2f} %\n")

    # # print the most common clipped sequences
    # sys.stderr.write("### most common 3' non-miRNA sequences detected\n")
    # for mirname, cdata in sorted(clipped.items()):
    #     # cdata = clipped[simple_name]
    #     sorted_counts = sorted(cdata.items(), key=lambda x: -x[1])
    #     for seq, n in sorted_counts[:50]:
    #         # if n > 9:
    #         sys.stderr.write(
    #             f"{'simple' if mirname in simple_names else 'full'}\t3'\t{mirname}\t{n}\t{100 * n/nclipped[mirname]:.2f} %\t{seq}\n"
    #         )
    #     if len(sorted_counts) >= 50:
    #         n_other = 0
    #         for seq, n in sorted_counts[50:]:
    #             n_other += n

    #         sys.stderr.write(
    #             f"{'simple' if mirname in simple_names else 'full'}\t3'\t{mirname}\t{n_other}\t{100 * n_other/nclipped[mirname]:.2f} %\tother\n"
    #         )

    # # print the most common clipped sequences
    # sys.stderr.write("### most common 5' non-miRNA sequences detected\n")
    # for mirname, cdata in sorted(clipped5.items()):
    #     # cdata = clipped[simple_name]
    #     for seq, n in sorted(cdata.items(), key=lambda x: -x[1])[:50]:
    #         if n > 9:
    #             sys.stderr.write(
    #                 f"{'simple' if mirname in simple_names else 'full'}\t5'\t{mirname}\t{n}\t{100 * n/nclipped5[mirname]:.2f} %\t{seq}\n"
    #             )
