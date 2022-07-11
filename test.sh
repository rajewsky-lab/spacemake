set -e

#rm project_df.csv > /dev/null

spacemake projects add_sample --project_id test \
    --sample_id sc_rnaseq_sample \
    --R1 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse   

spacemake projects add_sample --project_id test \
    --sample_id sc_rnaseq_sample_2 \
    --R1 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse  \
    --barcode_flavor visium

# with one bc file
spacemake projects add_sample --project_id test \
    --sample_id one_bc_file \
    --R1 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse  \
    --barcode_flavor visium \
    --puck visium

# with two bc files
spacemake projects add_sample --project_id test \
    --sample_id two_bc_files \
    --R1 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --R2 spacemake/data/test/visium_public_lane_joined_1m_R1.fastq.gz \
    --species mouse  \
    --barcode_flavor visium \
    --puck visium \
    --puck_barcode_file spacemake/data/test/test_bc1.csv spacemake/data/test/test_bc2.csv

# update sample
spacemake projects update_sample --project_id test \
    --sample_id two_bc_files \
    --investigator Test

spacemake projects merge_samples --merged_project_id test \
    --merged_sample_id test_merged \
    --project_id_list test \
    --sample_id_list one_bc_file two_bc_files

# this is expected to fail as has different barcode_flavor
spacemake projects merge_samples --merged_project_id test \
    --merged_sample_id test_merged_2 \
    --project_id_list test \
    --sample_id_list sc_rnaseq_sample two_bc_files

spacemake projects merge_samples --merged_project_id test \
    --merged_sample_id test_merged_2 \
    --project_id_list test \
    --sample_id_list sc_rnaseq_sample_2 two_bc_files
