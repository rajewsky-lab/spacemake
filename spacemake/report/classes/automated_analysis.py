import base64
import io
import logging
import os

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from jinja2 import Template
from scipy.sparse import csr_matrix
from spacemake.config import ConfigFile
from spacemake.project_df import ProjectDF
from spacemake.util import message_aggregation

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


logger_name = "spacemake.snakemake.scripts.n_intersect_sequences"
logger = logging.getLogger(logger_name)


def setup_parser(parser):
    """
    Set up command-line arguments for the script.

    :param parser: Argument parser object.
    :type parser: argparse.ArgumentParser
    :returns: Updated argument parser object.
    :rtype: argparse.ArgumentParser
    """
    # These allow to get the run_mode_variables.* from the config file
    # and pbf_metrics.px_by_um via project_df.puck_barcode_file_metrics
    parser.add_argument(
        "--project-id",
        type=str,
        help="the project_id in spacemake project_df.csv",
        required=True,
    )
    parser.add_argument(
        "--sample-id",
        type=str,
        help="the sample_id in spacemake project_df.csv",
        required=True,
    )
    parser.add_argument(
        "--run-mode",
        type=str,
        help="one of the run_mode from the project+sample in project_df.csv",
        required=True,
    )
    parser.add_argument(
        "--puck-barcode-file-id-qc",
        type=str,
        help="a path to the puck_barcode_file_id_qc",
        required=True,
    )
    # These have the paths to the input files used for generating plots
    parser.add_argument(
        "--obs-df",
        type=str,
        help="path to the csv file that stores properties per observation (per cell)",
        required=True,
    )
    parser.add_argument(
        "--var-df",
        type=str,
        help="path to the csv file that stores properties per variable (per gene)",
        required=True,
    )
    parser.add_argument(
        "--neighborhood-enrichment",
        type=str,
        help="path to the csv file with the neighborhood enrichment results",
        required=False,
        default="",
    )
    # These additionally configure the plotting functionality
    parser.add_argument(
        "--is-spatial",
        type=str,
        help="Whether the current sample is spatial or not",
        required=True,
        choices=["True", "False"],
        default="False",
    )
    # This specifies the UMI filter for the current report
    parser.add_argument(
        "--umi-cutoff",
        type=str,
        help="one of the UMI cutoffs from the --run-mode",
        required=True,
    )
    # This specifies where the output file will be generated
    parser.add_argument(
        "--out-html-report",
        type=str,
        help="where the HTML report will be saved",
        required=True,
    )

    return parser


def plot_metric(values, axis, nbins=100, color="#000000"):
    # decide linear or logarithmic scale
    min_difference = values.max() - values.min()

    hist, bins = np.histogram(values, bins=nbins)
    if np.abs(min_difference) < 100:
        axis.bar(bins[:-1], hist, color=color)

    else:
        logbins = np.logspace(np.log10(bins[0] + 1), np.log10(bins[-1]), nbins)
        axis.hist(values, bins=logbins, color=color)
        axis.set_xscale("log")

    axis.spines[["right", "top"]].set_visible(False)


def plot_to_base64(fig):
    my_stringIObytes = io.BytesIO()
    fig.savefig(my_stringIObytes, format="jpg")
    my_stringIObytes.seek(0)
    return base64.b64encode(my_stringIObytes.read()).decode()


def generate_automated_analysis_metadata(
    project_id: str,
    sample_id: str,
    run_mode: str,
    puck_barcode_file_id_qc: str,
    obs_df_path: str,
    var_df_path: str,
    neighborhood_enrichment: str,
    is_spatial: bool = False,
    umi_cutoff: int = 0,
):
    report = {
        "table_overview": [],
        "data_complete": False,
        "n_cells": 0,
        "umi_cutoff": 0,
        "n_genes": 0,
        "analyzed_clustering_resolutions": [],
        "is_spatial": False,
        "plots": [],
    }
    main_plots = {
        "histogram_metrics_spatial_units": "",
        "umi_spatial_unit_original": "",
        "umi_spatial_unit_scaled": "",
    }

    # Loading data
    obs_df = pd.read_csv(obs_df_path)
    var_df = pd.read_csv(var_df_path)

    # Loading project_df for metadata
    config = ConfigFile.from_yaml("config.yaml")
    project_df = ProjectDF("project_df.csv", config=config)

    # Create anndata for scanpy plotting
    empty_ad = csr_matrix((len(obs_df), len(var_df)), dtype=int)
    adata = ad.AnnData(empty_ad)
    adata.obs = obs_df
    adata.var = var_df
    adata.obsm["X_umap"] = adata.obs[["umap_0", "umap_1"]].values

    report["n_genes"] = len(var_df)
    report["n_cells"] = len(obs_df)

    report["umi_cutoff"] = umi_cutoff
    report["data_complete"] = True if np.any(obs_df.columns.str.startswith("leiden")) else False
    report["is_spatial"] = is_spatial

    report["plots"] = main_plots

    # Plot Histogram of metrics over spatial units
    fig, axes = plt.subplots(3, 2, figsize=(7, 6))
    plot_metric(obs_df["n_reads"], axes[0, 0], 100, clrs["reads"])
    axes[0, 0].set_ylabel("# of\nspatial units")
    axes[0, 0].set_xlabel("# of reads")
    plot_metric(obs_df["reads_per_counts"], axes[0, 1], 100, clrs["pcr"])
    axes[0, 1].set_ylabel("# of\nspatial units")
    axes[0, 1].set_xlabel("# of reads/UMI")
    plot_metric(obs_df["n_genes_by_counts"], axes[1, 0], 100, clrs["genes"])
    axes[1, 0].set_ylabel("# of\nspatial units")
    axes[1, 0].set_xlabel("# of genes")
    plot_metric(obs_df["total_counts"], axes[1, 1], 100, clrs["umis"])
    axes[1, 1].set_ylabel("# of\nspatial units")
    axes[1, 1].set_xlabel("# of UMIs")
    plot_metric(obs_df["pct_counts_mt"], axes[2, 0], 100, clrs["pct_counts_mt"])
    axes[2, 0].set_ylabel("# of\nspatial units")
    axes[2, 0].set_xlabel("% of mitochondrial counts")
    axes[2, 1].axis("off")
    plt.tight_layout()
    main_plots["histogram_metrics_spatial_units"] = plot_to_base64(fig)

    # Variables for spatial
    px_by_um, puck_bead_size = 1, 1
    x_breaks, y_breaks = None, None
    x_mm_breaks, y_mm_breaks = None, None
    puck_width_um = 1

    # Distribution of UMIs in 2D
    if is_spatial:
        # Loading metadata from spatial pucks and run mode
        puck_metrics = (
            project_df.get_puck_barcode_file_metrics(
                project_id=project_id,
                sample_id=sample_id,
                puck_barcode_file_id=puck_barcode_file_id_qc,
            ),
        )[0]

        puck_settings = project_df.get_puck_variables(
            project_id, sample_id, return_empty=True
        )

        run_mode_vars = project_df.config.get_run_mode(run_mode).variables

        print(puck_metrics)
        print(puck_settings)
        print(run_mode_vars)

        px_by_um = puck_metrics["px_by_um"]
        mesh_spot_diameter_um = run_mode_vars["mesh_spot_diameter_um"]
        meshed = run_mode_vars["mesh_data"]
        spot_diameter_um = puck_settings["spot_diameter_um"]

        adata.obsm["spatial"] = adata.obs[["x_pos", "y_pos"]].values

        # Set limits and axes units for the spatial plots
        x_limits = adata.obsm["spatial"][:, 0].min(), adata.obsm["spatial"][:, 0].max()
        y_limits = adata.obsm["spatial"][:, 1].min(), adata.obsm["spatial"][:, 1].max()
        puck_width_um = (x_limits[1] - x_limits[0]) / px_by_um

        ratio = (x_limits[1] - x_limits[0]) / (y_limits[1] - y_limits[0])

        scale_factor = 2 if puck_width_um < 3000 else 3
        mm_dist = max(10**scale_factor, round(puck_width_um / 3, scale_factor))
        mm_diff = mm_dist / 1000

        def_plot_bead_size = 0.5 if report["n_cells"] > 5000 else 0.75
        def_plot_bead_size = 0.1 if report["n_cells"] > 10000 else def_plot_bead_size
        def_plot_bead_size = 0.05 if report["n_cells"] > 25000 else def_plot_bead_size

        puck_bead_size = max(
            def_plot_bead_size, mesh_spot_diameter_um if meshed else spot_diameter_um
        )
        x_mm_breaks = np.arange(0, puck_width_um, mm_dist)
        x_mm_breaks = [f"{round(i, 1)} mm" for i in x_mm_breaks * mm_diff / mm_dist]
        y_mm_breaks = np.arange(0, puck_width_um / ratio, mm_dist)
        y_mm_breaks = [f"{round(i, 1)} mm" for i in y_mm_breaks * mm_diff / mm_dist]

        x_breaks = np.arange(x_limits[0], x_limits[1], px_by_um * mm_dist)
        y_breaks = np.arange(y_limits[0], y_limits[1], px_by_um * mm_dist)

        # Plot spatial UMI (original)
        fig, axes = plt.subplots(1, 1, figsize=(5, 5))
        sc.pl.spatial(
            adata,
            img_key=None,
            size=puck_bead_size * 1.5,
            spot_size=px_by_um,
            color="total_counts",
            ax=axes,
            show=False,
            cmap="inferno",
        )
        axes.spines[["right", "top"]].set_visible(False)
        axes.set_xticks(x_breaks)
        axes.set_xticklabels(x_mm_breaks)
        axes.set_yticks(y_breaks)
        axes.set_yticklabels(y_mm_breaks)
        axes.set_ylabel("")
        axes.set_xlabel("")

        main_plots["umi_spatial_unit_original"] = plot_to_base64(fig)

        # Plot spatial UMI (rescaled)

        fig, axes = plt.subplots(1, 1, figsize=(5, 5))
        sc.pl.spatial(
            adata,
            img_key=None,
            size=puck_bead_size * 1.5,
            spot_size=px_by_um,
            color="total_counts",
            ax=axes,
            show=False,
            cmap="inferno",
            vmax=np.quantile(adata.obs["total_counts"], 0.9),
        )
        axes.spines[["right", "top"]].set_visible(False)
        axes.set_xticks(x_breaks)
        axes.set_xticklabels(x_mm_breaks)
        axes.set_yticks(y_breaks)
        axes.set_yticklabels(y_mm_breaks)
        axes.set_ylabel("")
        axes.set_xlabel("")

        main_plots["umi_spatial_unit_scaled"] = plot_to_base64(fig)

    if report["data_complete"]:
        resolutions = np.unique(
            pd.Series(obs_df.columns[obs_df.columns.str.startswith("leiden")])
            .apply(lambda x: x.split("_")[1])
            .values
        )
        for i, resolution in enumerate(resolutions):
            _current_analyzed_clustering_resolution = {
        "resolution": "",
        "resolution_id": "",
        "plots": {"umap": "", "spatial": "", "neighbor_enrichment": ""},
    }

            _current_analyzed_clustering_resolution['resolution'] = f'{resolution} resolution'
            _current_analyzed_clustering_resolution['resolution_id'] = f'rid_{i}'

            # Plot UMAP for current resolution
            adata.obs[f"leiden_{resolution}"] = pd.Categorical(
                adata.obs[f"leiden_{resolution}"]
            )
            fig, axes = plt.subplots(1, 1, figsize=(5, 5))
            sc.pl.umap(adata, color=f"leiden_{resolution}", ax=axes, show=False)
            axes.spines[["right", "top", "left", "bottom"]].set_visible(False)
            axes.set_ylabel("UMAP 0")
            axes.set_xlabel("UMAP 1")
            _current_analyzed_clustering_resolution["plots"]["umap"] = plot_to_base64(
                fig
            )

            # Spatial plot and neighborhood enrichment (if spatial)
            if is_spatial:
                # Plot Spatial clusters
                fig, axes = plt.subplots(1, 1, figsize=(5, 5))
                sc.pl.spatial(
                    adata,
                    img_key=None,
                    size=puck_bead_size * 1.5,
                    spot_size=px_by_um,
                    color=f"leiden_{resolution}",
                    ax=axes,
                    show=False,
                )  # marker can be "o" or "h" if is hexagonal
                axes.spines[["right", "top"]].set_visible(False)
                axes.set_xticks(x_breaks)
                axes.set_xticklabels(x_mm_breaks)
                axes.set_yticks(y_breaks)
                axes.set_yticklabels(y_mm_breaks)
                axes.set_ylabel("")
                axes.set_xlabel("")
                _current_analyzed_clustering_resolution["plots"][
                    "spatial"
                ] = plot_to_base64(fig)

                # Plot Neighborhood enrichment
                nhood_enrich = pd.read_csv(neighborhood_enrichment)
                nhood_enrich_current_res = nhood_enrich[
                    nhood_enrich["resolution"].astype(str) == str(resolution)
                ]
                ne_data = np.zeros(
                    (
                        nhood_enrich_current_res["cluster_a"].max() + 1,
                        nhood_enrich_current_res["cluster_b"].max() + 1,
                    )
                )
                ne_data[
                    nhood_enrich_current_res["cluster_a"],
                    nhood_enrich_current_res["cluster_b"],
                ] = nhood_enrich_current_res["zscore"]

                fig, axes = plt.subplots(1, 1, figsize=(4, 4))
                plmat = axes.matshow(
                    ne_data, cmap="inferno", vmin=-50, vmax=100, origin="lower"
                )
                cbar = plt.colorbar(plmat, fraction=0.046)
                cbar.set_label("Neighbor enrichment score")
                axes.set_xlabel("cluster identity")
                axes.set_ylabel("cluster identity")
                axes.spines[["right", "top"]].set_visible(False)
                axes.xaxis.tick_bottom()
                plt.tight_layout()
                _current_analyzed_clustering_resolution["plots"][
                    "neighbor_enrichment"
                ] = plot_to_base64(fig)

            report["analyzed_clustering_resolutions"].append(
                _current_analyzed_clustering_resolution.copy()
            )

    report["table_overview"] = [
        {"name": "UMI filter", "value": umi_cutoff},
        {"name": "# genes in data", "value": report["n_genes"]},
        {"name": "# of spots in data", "value": report["n_cells"]},
        {"name": "median UMI", "value": obs_df["total_counts"].median()},
        {"name": "median genes", "value": obs_df["n_genes_by_counts"].median()},
        {"name": "puck width (um)", "value": puck_width_um},
    ]

    return report


def generate_html_report(data, template_file):
    with open(template_file, "r") as template_data:
        template_content = template_data.read()
        template = Template(template_content)

    html_report = template.render(report=data)

    return html_report


@message_aggregation(logger_name)
def cmdline():
    """cmdline."""
    import argparse

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="generate spacemake's 'automated analysis' with python",
    )
    parser = setup_parser(parser)

    args = parser.parse_args()
    template_file = os.path.join(absolute_path, "../templates/automated_analysis.html")

    report_metadata = generate_automated_analysis_metadata(
        args.project_id,
        args.sample_id,
        args.run_mode,
        args.puck_barcode_file_id_qc,
        args.obs_df,
        args.var_df,
        args.neighborhood_enrichment,
        args.is_spatial,
        args.umi_cutoff,
    )
    html_report = generate_html_report(report_metadata, template_file)

    with open(args.out_html_report, "w") as output:
        output.write(html_report)


if __name__ == "__main__":
    cmdline()