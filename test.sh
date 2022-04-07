set -e

#rm project_df.csv > /dev/null

spacemake projects add_sample --project_id test \
    --sample_id sc_rnaseq_sample \
    --R1 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse   

# with one bc file
spacemake projects add_sample --project_id test \
    --sample_id one_bc_file \
    --R1 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse  \
    --puck visium

# with two bc files
spacemake projects add_sample --project_id test \
    --sample_id two_bc_files \
    --R1 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse  \
    --puck visium \
    --puck_barcode_file spacemake/data/test/test_bc1.csv spacemake/data/test/test_bc2.csv
