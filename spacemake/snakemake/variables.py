# set root dir where the processed_data goes
project_dir = "projects/{project_id}"
config_path = "config.yaml"
project_df = "project_df.csv"

#################
# DIRECTORY STR #
#################
raw_data_root = project_dir + "/raw_data"
raw_data_illumina = raw_data_root + "/illumina"
raw_data_illumina_reads = raw_data_illumina + "/reads/raw"
raw_data_illumina_reads_reversed = raw_data_illumina + "/reads/bc_umi_tagged"
processed_data_root = project_dir + "/processed_data/{sample_id}"
processed_data_illumina = processed_data_root + "/illumina"

illumina_root = project_dir + "/processed_data/{sample_id}/illumina"
complete_data_root = illumina_root + "/complete_data"
data_root = illumina_root + "/{data_root_type}{downsampling_percentage}"
downsampled_data_prefix = illumina_root + "/downsampled_data"
downsampled_data_root = downsampled_data_prefix + "{downsampling_percentage}"

log_dir = complete_data_root + '/logs'
stats_dir = complete_data_root + '/stats'
plots_dir = complete_data_root + '/plots'

##############
# Demux vars #
##############
# Undetermined files pattern
# they are the output of bcl2fastq, and serve as an indicator to see if the demultiplexing has finished
demux_dir_pattern = "raw_data/demultiplex_data/{demux_dir}"
demux_indicator = demux_dir_pattern + "/indicator.log"

####################################
# FASTQ file linking and reversing #
####################################
reads_suffix = ".fastq.gz"

raw_reads_prefix = raw_data_illumina_reads + "/{sample_id}_R"
raw_reads_pattern = raw_reads_prefix + "{mate}" + reads_suffix
raw_reads_mate_1 = raw_reads_prefix + "1" + reads_suffix
raw_reads_mate_2 = raw_reads_prefix + "2" + reads_suffix

reverse_reads_prefix = raw_data_illumina_reads_reversed + "/{sample_id}_reversed_R"
reverse_reads_mate_1 = reverse_reads_prefix + "1" + reads_suffix

###############
# Fastqc vars #
###############
fastqc_root = raw_data_illumina + "/fastqc"
fastqc_pattern = fastqc_root + "/{sample_id}_R{mate}_fastqc.{ext}"
fastqc_ext = ["zip", "html"]

###
# splitting the reads
###

split_reads_root = complete_data_root + "/split_reads{polyA_adapter_trimmed}/"

split_reads_sam_names = [
    "plus_plus",
    "plus_minus",
    "minus_minus",
    "minus_plus",
    "plus_AMB",
    "minus_AMB",
]
split_reads_sam_pattern = split_reads_root + "{file_name}.sam"
split_reads_bam_pattern = split_reads_root + "{file_name}.bam"

split_reads_sam_files = [split_reads_root + x + ".sam" for x in split_reads_sam_names]

split_reads_strand_type = split_reads_root + "strand_type_num.txt"
split_reads_read_type = split_reads_root + "read_type_num.txt"

#######################
# post dropseq and QC #
#######################

qc_sheet = data_root + "/qc_sheets/qc_sheet_{sample_id}_{puck_barcode_file_id_qc}.html"
reads_type_out = split_reads_read_type
barcode_readcounts_suffix = "{polyA_adapter_trimmed}.txt.gz"
barcode_readcounts = complete_data_root + "/out_readcounts" + barcode_readcounts_suffix
barcode_readcounts_log = barcode_readcounts + ".log"
barcode_readcounts_prealigned = complete_data_root + "/out_readcounts_prealigned.txt.gz"
barcode_readcounts_prealigned_log = barcode_readcounts_prealigned + ".log"
strand_info = split_reads_strand_type

# united final.bam
top_barcodes_suffix = "{polyA_adapter_trimmed}.{n_beads}_beads.txt"
top_barcodes = complete_data_root + "/topBarcodes" + top_barcodes_suffix
top_barcodes_clean = complete_data_root + "/topBarcodesClean" + top_barcodes_suffix

spatial_barcodes = (
    complete_data_root
    + "/puck_barcode_files/spatialBarcodes_{puck_barcode_file_id}.txt"
)
parsed_spatial_barcodes = (
    complete_data_root
    + "/puck_barcode_files/spatial_barcodes_{puck_barcode_file_id}.csv"
)
parsed_spatial_barcodes_summary = (
    complete_data_root
    + "/puck_barcode_files/spatial_barcodes_summary_{puck_barcode_file_id}.csv"
)
parsed_spatial_barcodes_pc = (
    complete_data_root
    + "/puck_barcode_files/spatial_barcodes_puck_collection.csv"
)
stats_prealigned_spatial_barcodes = (
    complete_data_root
    + "/puck_barcode_files/stats_prealigned_spatial_barcodes_{puck_barcode_file_id}.csv"
)
puck_barcode_files_summary = complete_data_root + "/puck_barcode_files_summary.csv"
puck_count_barcode_matches_summary = complete_data_root + "/puck_count_barcode_matches.csv"
puck_count_prealigned_barcode_matches_summary = complete_data_root + "/puck_count_prealigned_barcode_matches.csv"

# dge creation
dge_root = data_root + "/dge"
dge_out_prefix = dge_root + "/dge"
dge_out_suffix = "{dge_type}{dge_cleaned}{polyA_adapter_trimmed}{mm_included}"
dge_out = (
    dge_out_prefix + dge_out_suffix + ".{n_beads}_beads_{puck_barcode_file_id}.txt.gz"
)
dge_out_summary = (
    dge_out_prefix
    + dge_out_suffix
    + ".{n_beads}_beads_{puck_barcode_file_id}.summary.txt"
)

# processed dge
h5ad_dge_suffix = "{is_external}.h5ad"
h5ad_dge_obs_suffix = "{is_external}.obs.csv"
dge_out_h5ad = (
    dge_out_prefix
    + dge_out_suffix
    + ".{n_beads}_beads_{puck_barcode_file_id}"
    + h5ad_dge_suffix
)
dge_out_h5ad_obs = (
    dge_out_prefix
    + dge_out_suffix
    + ".{n_beads}_beads_{puck_barcode_file_id}"
    + h5ad_dge_obs_suffix
)

# spatial dge
dge_spatial = (
    dge_out_prefix
    + dge_out_suffix
    + ".spatial_beads_{puck_barcode_file_id}"
    + h5ad_dge_suffix
)
dge_spatial_obs = (
    dge_out_prefix
    + dge_out_suffix
    + ".spatial_beads_{puck_barcode_file_id}"
    + h5ad_dge_obs_suffix
)

# spatial + collection dge
dge_spatial_collection = (
    dge_out_prefix
    + dge_out_suffix
    + ".spatial_beads_puck_collection"
    + h5ad_dge_suffix
)
dge_spatial_collection_obs = (
    dge_out_prefix
    + dge_out_suffix
    + ".spatial_beads_puck_collection"
    + h5ad_dge_obs_suffix
)

# spatial + meshed dge
dge_spatial_mesh_suffix = (
    ".spatial_beads.mesh_{spot_diameter_um}_{spot_distance_um}_{puck_barcode_file_id}"
)
dge_spatial_mesh_prefix = dge_out_prefix + dge_out_suffix + dge_spatial_mesh_suffix
dge_spatial_mesh = dge_spatial_mesh_prefix + h5ad_dge_suffix
dge_spatial_mesh_obs = dge_spatial_mesh_prefix + h5ad_dge_obs_suffix

# spatial + collection + meshed dge
dge_spatial_collection_mesh_suffix = (
    ".spatial_beads.mesh_{spot_diameter_um}_{spot_distance_um}_puck_collection"
)
dge_spatial_collection_mesh_prefix = dge_out_prefix + dge_out_suffix + dge_spatial_collection_mesh_suffix
dge_spatial_collection_mesh = dge_spatial_collection_mesh_prefix + h5ad_dge_suffix
dge_spatial_collection_mesh_obs = dge_spatial_collection_mesh_prefix + h5ad_dge_obs_suffix

dge_types = [
    ".exon",
    ".intron",
    ".all",
    ".Reads_exon",
    ".Reads_intron",
    ".Reads_all",
    "",
]

# kmer stats per position
kmer_stats_file = complete_data_root + "/kmer_stats/{kmer_len}mer_counts.csv"

# map paired-end to check errors
paired_end_prefix = complete_data_root + "/mapped_paired_end/"
paired_end_sam = paired_end_prefix + "{sample_id}_paired_end.sam"
paired_end_bam = paired_end_prefix + "{sample_id}_paired_end.bam"
paired_end_flagstat = paired_end_prefix + "{sample_id}_paired_end_flagstat.txt"
paired_end_log = paired_end_prefix + "{sample_id}_paired_end.log"
paired_end_mapping_stats = (
    paired_end_prefix + "{sample_id}_paired_end_mapping_stats.txt"
)

# automated analysis
automated_analysis_root = (
    data_root + "/automated_analysis/{run_mode}/umi_cutoff_{umi_cutoff}"
)
automated_report_prefix = (
    automated_analysis_root + "/{sample_id}_{puck_barcode_file_id_qc}_"
)
automated_report = automated_report_prefix + "automated_report.html"
automated_analysis_result_file = automated_report_prefix + "results.h5ad"

automated_analysis_processed_data_files = {
    "cluster_markers": "top10_cluster_markers.csv",
    "nhood_enrichment": "nhood_enrichment.csv",
    "obs_df": "obs_df.csv",
    "var_df": "var_df.csv",
}

# prepend automated_result_root
automated_analysis_processed_data_files = {
    key: automated_report_prefix + value
    for key, value in automated_analysis_processed_data_files.items()
}

# novosparc
novosparc_root = automated_analysis_root + "/novosparc"
novosparc_wildcards = (
    "/with_reference_rp_{reference_project_id}_"
    + "rs_{reference_sample_id}_rr_{reference_run_mode}_ru_{reference_umi_cutoff}"
)
novosparc_denovo_h5ad = novosparc_root + "/denovo.h5ad"
novosparc_with_reference_h5ad = novosparc_root + novosparc_wildcards + ".h5ad"

# in silico repo depletion
## This is the stderr output of bowtie2 after aligning to rRNA
ribo_depletion_log = complete_data_root + "/ribo_depletion_log.txt"
## This is the output of scripts/parse_ribo_log.py. Two lines:
# aligned_reads\t{n}
# input_reads\t{n}
parsed_ribo_depletion_log = complete_data_root + "/parsed_ribo_depletion_log.txt"

# #########################
#  dropseq rules and vars #
# #########################
tagged_bam = complete_data_root + "/unaligned_bc_tagged.bam"
tagged_bam_log = tagged_bam + ".log"
unassigned = complete_data_root + "/unaligned_bc_unassigned.bam"

# trim smart adapter from the reads
tagged_trimmed_bam = complete_data_root + "/unaligned_bc_tagged_trimmed.bam"

# trim polyA overheang if exists
tagged_polyA_adapter_trimmed_bam = (
    complete_data_root + "/unaligned_bc_tagged.polyA_adapter_trimmed.bam"
)

tagged_bam_pattern = (
    complete_data_root + "/unaligned_bc_tagged{polyA_adapter_trimmed}.bam"
)

# mapped reads
star_prefix = complete_data_root + "/star.{ref_name}."
star_log_file = complete_data_root + "/star.Log.final.out"
star_target_log_file = star_prefix + "Log.final.out"
star_tmp_dir = star_prefix + "tmp"

# final bam file
final_bam_suffix = "/final{polyA_adapter_trimmed}"
final_bam = complete_data_root + final_bam_suffix + ".bam"
bam_mm_included_pipe_suffix = "{dge_type}{dge_cleaned}{polyA_adapter_trimmed}.mm_included_{puck_barcode_file_id}.bam"
final_bam_mm_included_pipe = complete_data_root + "/final" + bam_mm_included_pipe_suffix

# downsampled bam
downsampled_bam_mm_included_pipe_suffix = "{dge_type}{dge_cleaned}{polyA_adapter_trimmed}.mm_included.bam"
downsampled_bam = (
    downsampled_data_root + "/final_downsampled{polyA_adapter_trimmed}.bam"
)
downsampled_bam_mm_included_pipe = (
    downsampled_data_root + "/final_downsampled" + downsampled_bam_mm_included_pipe_suffix
)
downsample_saturation_analysis = (
    downsampled_data_prefix + "/{project_id}_{sample_id}_{puck_barcode_file_id}_saturation_analysis.html"
)


# bt2_rRNA_index_basename = bt2_rRNA_index_dir + '/{species}_rRNA'
