import logging
import os

import numpy as np
import pandas as pd

# we use this instead of anndata, so we can load faster...
# safe assumption: the h5ad file is properly formatted...
import h5py

absolute_path = os.path.dirname(__file__)

cpalette = {
    "grey": "#999999",
    "light_orange": "#E69F00",
    "light_blue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "orange": "#D55E00",
    "pink": "#CC79A7",
}

clrs = {
    "umis": cpalette["light_orange"],
    "genes": cpalette["light_blue"],
    "reads": cpalette["green"],
    "pcr": cpalette["pink"],
    "pct_counts_mt": "black",
}

nucl_clrs = {
    "A": "#F5C900",
    "C": "#F55D59",
    "T": "#3AA861",
    "G": "#7772F5",
    "N": "#999999",
}

SAMPLEINFO_VARS = [
    "species",
    "sequencing_date",
    "investigator",
    "experiment",
]
SPATIAL_METRICS = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "n_reads",
    "reads_per_counts",
    "n_joined",
    "exact_entropy",
    "exact_compression",
]
SPATIAL_METRICS_TITLES = {
    "n_genes_by_counts": "# of genes per spatial unit",
    "total_counts": "# of UMIs per spatial unit",
    "pct_counts_mt": "# % mt counts per spatial unit",
    "n_reads": "# of reads per spatial unit",
    "reads_per_counts": "reads/UMI per spatial unit",
    "n_joined": "# beads joined per spatial unit",
    "exact_entropy": "Shannon entropy per spatial unit",
    "exact_compression": "barcode length after compression per spatial unit",
}


STRTOBOOL = {"False": False, "True": True}


logger_name = "spacemake.report.qc_sequencing"
logger = logging.getLogger(logger_name)


def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename, "rb") as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size)).decode(encoding="utf-8")
            remaining_size -= buf_size
            lines = buffer.split("\n")
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first
                if buffer[-1] != "\n":
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment


def read_star_log_file(log_file):
    file = os.path.join(log_file)

    if not os.path.isfile(file):
        return (1, 1)

    star_stats = {
        "input_reads": 0,
        "uniq_mapped_reads": 0,
        "multi_mapped_reads": 0,
        "too_many_mapped_reads": 0,
        "too_many_mapped_reads": 0,
        "unmapped_too_short": 0,
        "unmapped_other": 0,
        "chimeric": 0,
    }

    log_name_stat = {
        "Number of input reads": "input_reads",
        "Uniquely mapped reads number": "uniq_mapped_reads",
        "Number of reads mapped to multiple loci": "multi_mapped_reads",
        "Number of reads mapped to too many loci": "too_many_mapped_reads",
        "Number of reads unmapped: too many mismatches": "unmapped_mismatches",
        "Number of reads unmapped: too short": "unmapped_too_short",
        "Number of reads unmapped: other": "unmapped_other",
        "Number of chimeric reads": "chimeric",
        "Average mapped length": "avg_mapped_length",
    }

    with open(file) as fi:
        idx = 0
        for line in fi:
            _log_id = line.strip().split("|")[0].strip()
            if _log_id in log_name_stat.keys():
                val_str = line.strip().split("|")[1].strip()
                key = log_name_stat[_log_id]
                value = float(val_str) if val_str != "NA" else np.nan
                if _log_id == "Average mapped length":
                    star_stats[key] = value
                else:
                    star_stats[key] = round(value / 1e6, 2)
                # Do not convert to millions the Average mapped length

            idx += 1

    return star_stats


def read_split_reads_file(split_reads_file):
    file = os.path.join(split_reads_file)

    read_types = dict()
    with open(file) as fi:
        for line in fi:
            line = line.strip().split(" ")
            read_types[line[0]] = round(int(line[1]) / 1e6, 2)

    return read_types


def read_bowtie_log_file(log_file):
    file = os.path.join(log_file)
    bowtie_stats = {
        "input_reads": 0,
        "unique_aligned": 0,
        "multimapper": 0,
        "unaligned": 0,
    }
    log_line_stat = {
        5: "input_reads",
        3: "unaligned",
        2: "unique_aligned",
        1: "multimapper",
    }
    max_line = max(log_line_stat.keys())

    idx = 0
    for line in reverse_readline(file):
        if idx in log_line_stat.keys():
            bowtie_stats[log_line_stat[idx]] = round(
                int(line.strip().split(" ")[0]) / 1e6, 2
            )
        idx += 1
        if max_line < idx:
            break

    return bowtie_stats


def generate_table_mapping_statistics(
    complete_data_root: str, split_reads_read_type: str
):
    # Initialize empty lists to store the filenames
    bowtie_log_files = []
    star_log_files = []

    # Iterate over the files in the folder
    for filename in os.listdir(complete_data_root):
        # Check if the file ends with .bam.log
        if filename.endswith(".cram.log") and "final" not in filename:
            bowtie_log_files.append(filename)
        # Check if the file ends with .Log.final.out
        elif filename.endswith(".Log.final.out") and filename != "star.Log.final.out":
            star_log_files.append(filename)

    bowtie_logs = []
    for bowtie_log_file in bowtie_log_files:
        bowtie_log = read_bowtie_log_file(
            os.path.join(complete_data_root, bowtie_log_file)
        )
        bowtie_log["name"] = bowtie_log_file.split(".")[0]
        bowtie_log["mapper"] = "bowtie2"
        bowtie_logs.append(bowtie_log)

    star_logs = []
    for star_log_file in star_log_files:
        star_log = read_star_log_file(os.path.join(complete_data_root, star_log_file))
        star_log["name"] = star_log_file.split(".")[1]
        star_log["mapper"] = "STAR"
        star_logs.append(star_log)

    # we sort the mapping statistics by the number of input reads, merge into a single list
    all_logs = star_logs + bowtie_logs

    idx_log = np.argsort([log["input_reads"] for log in all_logs])[::-1]
    all_logs = [all_logs[idx] for idx in idx_log]

    reads_type = read_split_reads_file(split_reads_read_type)
    reads_type["name"] = "reads_type"
    all_logs += [reads_type]

    return all_logs


def create_run_modes_df(run_modes, project_df):
    all_variables = {}

    # Collect variables from all run modes
    for run_mode in run_modes:
        run_mode_vars = project_df.config.get_run_mode(run_mode)
        for var_name, var_value in run_mode_vars.variables.items():
            if var_name not in all_variables:
                all_variables[var_name] = {}
            all_variables[var_name][run_mode] = var_value

    run_modes_df = pd.DataFrame.from_dict(all_variables, orient="index")

    for mode in run_modes:
        if mode not in run_modes_df.columns:
            run_modes_df[mode] = None

    return run_modes_df


def create_mapping_stats_df(split_reads_read_type, data_root):
    if isinstance(split_reads_read_type, str):
        split_reads_read_type = [split_reads_read_type]

    all_stats = []

    for s in split_reads_read_type:
        strategy_name = os.path.basename(os.path.dirname(s))

        stats = generate_table_mapping_statistics(data_root, s)

        if not isinstance(stats, pd.DataFrame):
            stats = pd.DataFrame(stats)

        stats["Mapping Strategy"] = strategy_name

        all_stats.append(stats)

    combined_stats = pd.concat(all_stats, axis=0)
    print(combined_stats)

    # Extract only the rows for rRNA and STAR mapper (genome)
    rRNA_data = combined_stats[combined_stats["name"] == "rRNA"]
    STAR_data = combined_stats[combined_stats["mapper"] == "STAR"]
    reads_type_data = combined_stats[combined_stats["name"] == "reads_type"]

    # Get the mapping strategy (assuming it's the same for all)
    mapping_strategy = combined_stats["Mapping Strategy"].iloc[0]

    # Create the DataFrame structure with index as metrics and a single column for mapping_strategy
    metrics = [
        "input_reads",
        "uniq_mapped_reads",
        "as.cds",
        "as.utr",
        "intronic",
        "intergenic",
        "ambiguous",
        "avg_mapped_length",
        "multi_mapped_reads",
        "unmapped_too_short",
        "mapped_to_rRNA",
    ]

    # Create empty DataFrame with metrics as index and mapping_strategy as column
    result_df = pd.DataFrame(index=metrics, columns=[mapping_strategy])

    # Variables to store
    input_reads = 0

    # Get values from STAR mapper
    if not STAR_data.empty:
        star_row = STAR_data.iloc[0]

        # Extract the metrics we want
        input_reads = star_row.get("input_reads", 0)
        uniq_mapped_reads = star_row.get("uniq_mapped_reads", 0)
        multi_mapped_reads = star_row.get("multi_mapped_reads", 0)
        unmapped_too_short = star_row.get("unmapped_too_short", 0)
        avg_mapped_length = star_row.get("avg_mapped_length", 0)

        # Calculate percentages
        if input_reads > 0:
            uniq_pct = (uniq_mapped_reads / input_reads) * 100
            multi_pct = (multi_mapped_reads / input_reads) * 100
            too_short_pct = (unmapped_too_short / input_reads) * 100
        else:
            uniq_pct = multi_pct = too_short_pct = 0

        # Add to result DataFrame
        result_df.loc["input_reads", mapping_strategy] = f"{input_reads:.2f}"
        result_df.loc["uniq_mapped_reads", mapping_strategy] = (
            f"{uniq_mapped_reads:.2f} ({uniq_pct:.1f}%)"
        )
        result_df.loc["avg_mapped_length", mapping_strategy] = (
            f"{avg_mapped_length:.2f}"
        )
        result_df.loc["multi_mapped_reads", mapping_strategy] = (
            f"{multi_mapped_reads:.2f} ({multi_pct:.1f}%)"
        )
        result_df.loc["unmapped_too_short", mapping_strategy] = (
            f"{unmapped_too_short:.2f} ({too_short_pct:.1f}%)"
        )

    # Get values from reads_type
    if not reads_type_data.empty:
        type_row = reads_type_data.iloc[0]

        # Extract the metrics
        coding = type_row.get("CODING", 0)
        utr = type_row.get("UTR", 0)
        intronic = type_row.get("INTRONIC", 0)
        intergenic = type_row.get("INTERGENIC", 0)
        ambiguous = type_row.get("AMB", 0)

        # Calculate percentages
        if input_reads > 0:
            coding_pct = (coding / input_reads) * 100
            utr_pct = (utr / input_reads) * 100
            intronic_pct = (intronic / input_reads) * 100
            intergenic_pct = (intergenic / input_reads) * 100
            ambiguous_pct = (ambiguous / input_reads) * 100
        else:
            coding_pct = utr_pct = intronic_pct = intergenic_pct = ambiguous_pct = 0

        # Add to result DataFrame
        result_df.loc["as.cds", mapping_strategy] = f"{coding:.2f} ({coding_pct:.1f}%)"
        result_df.loc["as.utr", mapping_strategy] = f"{utr:.2f} ({utr_pct:.1f}%)"
        result_df.loc["intronic", mapping_strategy] = (
            f"{intronic:.2f} ({intronic_pct:.1f}%)"
        )
        result_df.loc["intergenic", mapping_strategy] = (
            f"{intergenic:.2f} ({intergenic_pct:.1f}%)"
        )
        result_df.loc["ambiguous", mapping_strategy] = (
            f"{ambiguous:.2f} ({ambiguous_pct:.1f}%)"
        )

    # Get values from rRNA
    if not rRNA_data.empty:
        rRNA_row = rRNA_data.iloc[0]
        rRNA_reads = rRNA_row.get("unique_aligned", 0)

        # Calculate percentage
        if input_reads > 0:
            rRNA_pct = (rRNA_reads / input_reads) * 100
        else:
            rRNA_pct = 0

        # Add to result DataFrame
        result_df.loc["mapped_to_rRNA", mapping_strategy] = (
            f"{rRNA_reads:.2f} ({rRNA_pct:.1f}%)"
        )

    # Reset the index to make it a regular DataFrame column
    result_df = result_df.rename_axis("Metric").reset_index()

    return result_df


def create_summary_beads_df(run_modes_adatas):
    all_variables = {
        "median_genes": {},
        "median_pcr": {},
        "median_reads": {},
        "median_umis": {},
        "n_beads": {},
        "sum_reads": {},
    }

    # Collect variables from all run modes
    for run_mode, adata in run_modes_adatas.items():
        with h5py.File(adata) as _adata_f:
            all_variables["median_genes"][run_mode] = np.median(
                _adata_f["obs/n_genes_by_counts"]
            )
            all_variables["median_pcr"][run_mode] = np.median(
                _adata_f["obs/reads_per_counts"]
            )
            all_variables["median_reads"][run_mode] = np.median(_adata_f["obs/n_reads"])
            all_variables["median_umis"][run_mode] = np.median(
                _adata_f["obs/total_counts"]
            )
            all_variables["n_beads"][run_mode] = len(_adata_f["obs/total_counts"])
            all_variables["sum_reads"][run_mode] = np.round(
                _adata_f["obs/n_reads"][:].sum() / 1e6, 3
            )

    run_modes_df = pd.DataFrame.from_dict(all_variables, orient="index")

    for mode in run_modes_adatas.keys():
        if mode not in run_modes_df.columns:
            run_modes_df[mode] = None

    return run_modes_df


def create_sample_info_df(project_df, project_id, sample_id, puck_barcode_file_id_qc):
    sample_info = project_df.get_sample_info(project_id, sample_id)
    sample_info_df = pd.DataFrame(
        {
            "project_id": [project_id],
            "sample_id": [sample_id],
            "puck_barcode_file_id": [puck_barcode_file_id_qc],
            "species": [sample_info["species"]],
            "sequencing_date": [sample_info["sequencing_date"]],
            "investigator": [sample_info["investigator"]],
            "experiment": [sample_info["experiment"]],
        }
    )

    return sample_info_df


def create_metrics_table_df(adata, umi_cutoff):
    metrics_dict = {
        "UMI filter": [umi_cutoff],
        "Number of genes in data": [len(adata.var_names)],
        "Number of spots in data": [len(adata)],
        "Median UMI": [np.median(adata.obs["total_counts"])],
        "Median Genes": [np.median(adata.obs["n_genes_by_counts"])],
    }

    if "puck_variables" in adata.uns and "width_um" in adata.uns["puck_variables"]:
        metrics_dict["Puck width (µm)"] = [adata.uns["puck_variables"]["width_um"]]

    metrics_table_df = pd.DataFrame(metrics_dict)

    return metrics_table_df


def generate_deciled_data(values):
    # Group by '_downsample_pct_report' and perform the necessary calculations
    def calculate_deciles(group):
        group["cumsum_reads"] = group["n_reads"].cumsum()
        group["decile_limit"] = group["n_reads"].sum() / 10
        group["decile"] = (group["cumsum_reads"] / group["decile_limit"]).floordiv(
            1
        ) + 1
        return group.loc[group["decile"] < 11]

    # Group by '_downsample_pct_report' and apply the calculate_deciles function
    decile_dat = (
        values.groupby("_downsample_pct_report")
        .apply(calculate_deciles)
        .reset_index(drop=True)
    )

    # Group by 'percentage' and 'decile' and calculate medians and counts
    decile_dat = (
        decile_dat.groupby(["_downsample_pct_report", "decile"])
        .agg(
            {
                "n_reads": "median",
                "n_genes_by_counts": "median",
                "reads_per_counts": "median",
                "total_counts": "median",
                "cell_bc": "count",
            }
        )
        .reset_index()
    )

    # Melt the DataFrame to long format
    decile_dat = pd.melt(
        decile_dat,
        id_vars=["_downsample_pct_report", "decile"],
        var_name="observation",
        value_name="value",
    )

    # Convert 'decile' and '_downsample_pct_report' to appropriate data types
    decile_dat["decile"] = decile_dat["decile"].astype("category")
    decile_dat["_downsample_pct_report"] = decile_dat["_downsample_pct_report"].astype(
        int
    )

    mapping_dict = {
        "n_reads": "median_reads",
        "n_genes_by_counts": "median_genes",
        "reads_per_counts": "median_pcr",
        "total_counts": "median_umis",
        "cell_bc": "n_beads",
    }
    decile_dat["observation"] = decile_dat["observation"].replace(mapping_dict)

    return decile_dat


def load_dge_summary_downsampling(dge_summaries, run_mode, downsample_pcts, puck_barcode_file_id):
    obs_df = pd.DataFrame()
    for downsample_pct in downsample_pcts:
        _path = dge_summaries[f"downsampled_dge_summary.{run_mode}.{downsample_pct}.{puck_barcode_file_id}"]
        if isinstance(_path, list):
            _path = _path[0]
        _obs_df = pd.read_csv(_path)
        _obs_df["_downsample_pct_report"] = downsample_pct
        obs_df = pd.concat([obs_df, _obs_df])

    return obs_df