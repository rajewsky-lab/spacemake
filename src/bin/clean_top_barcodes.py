import pysam
import difflib

# we need to reverse it
optical_primer = 'GAATCACGATACGTACACCA'[::-1]
optical_primer_len = len(optical_primer)

nucl_stretches = ['TTTTTT', 'AAAAAAAA', 'CCCCCCCC', 'GGGGGGGG']

with open(snakemake.input[0], 'r') as fi, open(snakemake.output[0], 'w') as fo:
    for barcode in fi:
        barcode = barcode.strip()
        barcode_len = len(barcode)
        
        # clean up TAG=XC artifact
        if barcode == 'TAG=XC':
            continue
        
        matcher = difflib.SequenceMatcher(None, optical_primer, barcode)
        
        pos_optical_primer, pos_barcode, kmer_len = matcher.find_longest_match(0, optical_primer_len, 0, barcode_len)
        
        # if overlap with barcode is bigger than 4, and the overlap is at the end, skip
        if kmer_len > 3 and pos_barcode + kmer_len == barcode_len:
            continue
        
        # if overlap at least 7, anywhere, skip
        if kmer_len > 6:
            continue
        
        # if any of the nucl stretches is in the barcode, skip
        if any([stretch in barcode for stretch in nucl_stretches]):
            continue
        
        # write line to file
        _ = fo.write(barcode + '\n')
