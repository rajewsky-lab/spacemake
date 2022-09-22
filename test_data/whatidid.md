# Creat reasonably sized test data for a number of SPACEMAKE low-level features

## Reference sequences

### Genome

I chose to extract chr22 (the smallest) from hg38 using some custom script.

### Contaminating rRNA

I chose to use my index of abundant human rRNA, tRNA, etc. sequences, hand built for other projects

### miRNA

I use miRGeneDB 2.1

## Annotation

### Genome

Extracting only chr22 from GENCODE v38

```    
    grep chr22 /data/rajewsky/annotation/GRCh38/gencode.v38.annotation.gtf > gencode.v38.chr22.gtf
```

### miRNA

Created a pseudo-annotation that assigns each mirgenedb entry to the miRNA gene name

# Build the species references

```
    spacemake init --dropseq_tools=/data/rajewsky/shared_bins/Drop-seq_tools-2.4.0

    spacemake config add_species \
        --name=test_hsa \
        --reference=genome \
        --sequence=hg38_chr22.fa \
        --annotation=gencode.v38.chr22.gtf
    
    spacemake config add_species \
        --name=test_hsa \
        --reference=rRNA \
        --sequence=rRNA_hsa.fa

    spacemake config add_species \
        --name=test_hsa \
        --reference=miRNA \
        --sequence=mirgenedb.hsa.mature.fa \
        --annotation=mirgenedb.hsa.mature.gtf
```

After this the `config.yaml` should contain the following section:

```
    species:
        test_hsa:
            rRNA:
                BT2_index: null
                STAR_index_dir: null
                annotation: ''
                sequence: rRNA_hsa.fa
            genome:
                BT2_index: null
                STAR_index_dir: null
                annotation: gencode.v38.chr22.gtf
                sequence: hg38_chr22.fa
            miRNA:
                BT2_index: null
                STAR_index_dir: null
                annotation: mirgenedb.hsa.mature.gtf
                sequence: mirgenedb.hsa.mature.fa
```

## Create some dummy sample with a handful of hand-picked reads

Hand-crafted test_reads.R1.fastq.gz and test_reads.R2.fastq.gz to contain rRNA, two miRNA reads, and a few mRNA (coding, alterntaive UTR), as well as an intergenic and a random sequence. The cell barcodes are identical and the UMI is simply counting up.

## Add the sample and have it run

```
    spacemake projects add_sample \
        --species=test_hsa \
        --project_id=test \
        --sample_id=test_01 \
        --R1=test_reads.R1.fastq.gz \
        --R2=test_reads.R2.fastq.gz \
        --map_strategy="bowtie2:aRNA->bowtie2:miRNA->STAR:genome:final"
``` 

Small test-data have been added to spacemake/test_data and the steps to add sequences,annotation and samples were added to spacemake/spacemake/unittests.py

# Part II: More extensive test-data to stress-test annotation code

Nikos has kindly extracted a couple dozen uniquely mapping reads in all sorts of relationship to a single gene model (UTR, coding, intronic, ...) from chr22.

```
    cat /data/rajewsky/home/nkarais/murphy/fc_sts/collect_reads_chr22/unittest_reads_2.fastq | gzip -c > reads_chr22_R2.fastq.gz

    cat /data/rajewsky/home/nkarais/murphy/fc_sts/collect_reads_chr22/unittest_reads_1.fastq | gzip -c > reads_chr22_R1.fastq.gz
```

We are going to extract the entire genic loci (plus some intergenic flank) that are required to map these reads, and do the same for the gencode annotation.

After that, we can have test-cases for the annotation code! (towards independence of the dropseq tools for SPACEMAKE 2.x).

## Extract the gene regions and corresponding annotation

```
    python gene_loci_from_gtf.py | sort -grk 4 > chr22_gene_bounds.csv
    python make_chr22_test_data.py 
```
