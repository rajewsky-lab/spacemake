import pytest
import spacemake.quant as quant
import os
test_bundles = [
    ("empty", 'ACGTACGT', 'AAAAA', [], ("-", set())),
    ("uniq", 'ACGTACGT', 'AAAAC', [('chr1', '+', 'A', 'N', 0)], ("A", {'counts', 'reads', 'exonic_counts', 'exonic_reads'})),
    ("mm_gene_intergenic", 'ACGTACGT', 'AAAAG', [
        ('chr1', '+', 'A', 'C', 0),
        ('chr2', '-', '-', '-', 0),
        ], ("A", {'counts', 'reads', 'exonic_counts', 'exonic_reads'})),
    ("mm_geneC_geneI", 'ACGTACGT', 'AAACA', [
        ('chr1', '+', 'A', 'C', 0),
        ('chr2', '-', 'B', 'I', 0),
        ], ("A", {'counts', 'reads', 'exonic_counts', 'exonic_reads'})),
    ("mm_complex", 'ACGTACGT', 'ACACA', [
        ('chr1', '+', 'A', 'C|N|I,I', 0),
        ('chr2', '-', 'B', 'n', 0),
        ], ("A", {'counts', 'reads', 'intronic_counts', 'intronic_reads'})),
    ("mm_complex2", 'ACGTACGT', 'ATACA', [
        ('chr1', '+', 'A', 'N,C|N|I,U', 0),
        ('chr2', '-', 'B', 'n', 0),
        ], ("A", {'counts', 'reads', 'exonic_counts', 'exonic_reads'})),
    ("mm_ambig", 'ACGTACGT', 'AAACG', [
        ('chr1', '+', 'A', 'C', 0),
        ('chr2', '-', 'B', 'C', 0),
        ], (None, set())),
    ("mm_geneI_geneN", 'ACGTACGT', 'AAACT', [
        ('chr1', '+', 'A,B', 'I,n', 0),
        ], ("A", {'counts', 'reads', 'intronic_counts', 'intronic_reads'})),
    ("uniq_dup", 'ACGTACGT', 'AAAAC', [('chr1', '+', 'A', 'N', 0)], ("A", {'reads', 'exonic_reads'})),
    ("uniq", 'NNNNNNNN', 'AAAAC', [('chr1', '+', 'A', 'N', 0)], ("A", {'counts', 'reads', 'exonic_counts', 'exonic_reads'})),
]

def test_default_counter():
    counter = quant.DefaultCounter()
    sm_dir = os.path.dirname(__file__)
    dge = quant.DGE(channels=counter.channels, cell_bc_allowlist=sm_dir + '/../test_data/simple_allowlist.txt')

    fail = []
    for name, cb, umi, bundle, expect in test_bundles:
        res = counter.process_bam_bundle(cell=cb, umi=umi, bundle=bundle)
        gene, channels = res
        print(res == expect, res, expect)
        if res != expect:
            fail.append((name, bundle, res, expect))
        
        dge.add_read(gene, cb, channels)
    
    channel_d, obs, var = dge.make_sparse_arrays()

    adata = quant.DGE.sparse_arrays_to_adata(channel_d, obs, var)
    print(adata)

