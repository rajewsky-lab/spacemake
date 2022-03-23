################################
# Final output file generation #
################################
def get_output_files(
    pattern,
    projects=[],
    samples=[],
    filter_merged=False,
    run_on_external=True,
    **kwargs,
):
    out_files = []
    df = project_df.df

    if projects != [] or samples != []:
        project_df.assert_index_value(projects, "project_id")

        project_df.assert_index_value(samples, "sample_id")

        ix = project_df.get_ix_from_project_sample_list(
            project_id_list=projects, sample_id_list=samples
        )

        df = df.loc[ix]

    if filter_merged:
        df = df.loc[~df.is_merged]

    for index, row in df.iterrows():
        is_external = project_df.is_external(
            project_id=index[0], sample_id=index[1]
        )

        has_dge = project_df.has_dge(
            project_id=index[0], sample_id=index[1]
        )

        if not has_dge:
            # if there is no dge, skip
            continue
        
        if not run_on_external and is_external:
            # if we do not want to run this on external samples
            # and sample is external, skip
            continue

        for run_mode in row["run_mode"]:
            run_mode_variables = project_df.config.get_run_mode(run_mode).variables

            out_files = out_files + expand(
                pattern,
                project_id=index[0],
                sample_id=index[1],
                puck_id=row["puck_id"],
                run_mode=run_mode,
                umi_cutoff=run_mode_variables["umi_cutoff"],
                **kwargs,
            )

    return out_files

def get_all_dges(wildcards):
    df = project_df.df

    dges = []

    for index, row in df.iterrows():
        for run_mode in row["run_mode"]:
            if project_df.has_dge(project_id=index[0], sample_id=index[1]):
                dges.append(
                    get_dge_from_run_mode(
                        project_id=index[0],
                        sample_id=index[1],
                        run_mode=run_mode,
                        data_root_type="complete_data",
                        downsampling_percentage="",
                    )["dge"],
                )

    return dges


def get_raw_dge(wildcards):
    is_external = project_df.is_external(
        project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )

    has_dge = project_df.has_dge(
        project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )

    if has_dge:
        out_files = {}
    else:
        return []

    if is_external:
        out_files["dge"] = project_df.get_metadata(
            "dge", project_id=wildcards.project_id, sample_id=wildcards.sample_id
        )
    else:
        out_files["dge"] = dge_out
        out_files["dge_summary"] = dge_out_summary

    return out_files


def get_reads(wildcards):
    ###
    # R1 and R2 for demultiplexed reads will return none
    ###
    reads = project_df.get_metadata(
        "R" + wildcards.mate, sample_id=wildcards.sample_id, project_id=wildcards.project_id
    )
    if reads is None or reads == []:
        return ["none"]
    else:
        # reads already
        return reads


# barcode flavor parsing and query functions
class dotdict(dict):
    """dot.notation access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def parse_barcode_flavors(
    config,
    bc_default_settings=dict(
        bc1_ref="",
        bc2_ref="",
        cell_raw="None",
        score_threshold=0.0,
        min_opseq_score=22,
        bam_tags="CR:{cell},MI:{UMI}",
    ),
):
    """
    Reads the 'barcode_flavor' from 'knowledge' section of the config.yaml.
    parses and gathers the settings for barcode flavors
    """
    preprocess_settings = {}
    for flavor, flavor_settings in config["barcode_flavors"].items():
        # for each flavor, also retrieve the configuration
        # first make a copy of the default values
        d = dict(bc_default_settings)
        d.update(flavor_settings)
        preprocess_settings[flavor] = dotdict(d)

    res = dotdict(
        dict(
            preprocess_settings=preprocess_settings,
        )
    )

    return res


# all barcode flavor info from config.yaml
# is kept here for convenient lookup
bc_flavor_data = parse_barcode_flavors(config)


def get_bc_preprocess_settings(wildcards):
    """
    This function will return a dictionary of information
    on the read1 preprocessing, according to barcode_flavor
    """
    flavor = project_df.get_metadata(
        "barcode_flavor", project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )
    if flavor not in bc_flavor_data.preprocess_settings:
        raise Exception(flavor)

    settings = bc_flavor_data.preprocess_settings[flavor]

    return settings


def get_demux_indicator(wildcards):
    demux_dir = project_df.get_metadata(
        "demux_dir", sample_id=wildcards.sample_id, project_id=wildcards.project_id
    )

    return expand(demux_indicator, demux_dir=demux_dir)


def get_star_input_bam(wildcards):
    if wildcards.polyA_adapter_trimmed == ".polyA_adapter_trimmed":
        return {"reads": tagged_polyA_adapter_trimmed_bam}
    else:
        return {"reads": tagged_bam}


def get_final_bam(wildcards):
    is_merged = project_df.get_metadata(
        "is_merged", project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )

    if is_merged:
        return [final_merged_bam]
    else:
        return [final_bam]


def get_dge_input_bam(wildcards):
    if wildcards.data_root_type == "complete_data":
        final_bam_pipe = final_bam_mm_included_pipe
        final_bam = get_final_bam(wildcards)
    elif wildcards.data_root_type == "downsampled_data":
        final_bam_pipe = downsampled_bam_mm_included_pipe
        final_bam = downsampled_bam

    if wildcards.mm_included == ".mm_included":
        out = {"reads": final_bam_pipe}
    else:
        out = {"reads": final_bam}

    return out

def get_species_genome(wildcards):
    # This function will return the genome of a sample
    #    - annotation (.gtf file)
    #    - genome (.fa file)
    if "species" not in wildcards.keys():
        species = project_df.get_metadata(
            "species", project_id=wildcards.project_id, sample_id=wildcards.sample_id
        )
    else:
        species = wildcards.species

    files = project_df.config.get_variable("species", name=species)

    return [files['genome']]

def get_species_annotation(wildcards):
    # This function will return the genome of a sample
    #    - annotation (.gtf file)
    #    - genome (.fa file)
    if "species" not in wildcards.keys():
        species = project_df.get_metadata(
            "species", project_id=wildcards.project_id, sample_id=wildcards.sample_id
        )
    else:
        species = wildcards.species

    files = project_df.config.get_variable("species", name=species)

    return [files['annotation']]

def get_species_genome_annotation(wildcards):
    # This function will return 2 things required by STAR:
    #    - annotation (.gtf file)
    #    - genome (.fa file)
    if "species" not in wildcards.keys():
        species = project_df.get_metadata(
            "species", project_id=wildcards.project_id, sample_id=wildcards.sample_id
        )
    else:
        species = wildcards.species

    return {'genome': species_genome.format(species=species),
            'annotation': species_annotation.format(species=species)}


def get_star_index(wildcards):
    # This function will return 1 things required by STAR:
    #    - index directory
    species = project_df.get_metadata(
        "species", project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )
    species_data = project_df.config.get_variable("species", name=species)

    if 'STAR_index_dir' in species_data:
        return {'index': species_data['STAR_index_dir']}
    else:
        return {"index": expand(star_index, species=species)[0]}

def get_rRNA_genome(wildcards):
    files = project_df.config.get_variable("species", name=wildcards.species)

    return [files["rRNA_genome"]]


def get_bt2_rRNA_index(wildcards):
    species = project_df.get_metadata(
        "species", project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )

    files = project_df.config.get_variable("species", name=species)

    if "rRNA_genome" in files:
        return {"index": expand(bt2_rRNA_index_dir, species=species)[0]}

    return []


def get_run_modes_from_sample(project_id, sample_id):
    run_mode_names = project_df.get_metadata(
        "run_mode", project_id=project_id, sample_id=sample_id
    )

    run_modes = {}

    for run_mode in run_mode_names:
        run_modes[run_mode] = project_df.config.get_run_mode(run_mode).variables

    return run_modes


def get_dge_extra_params(wildcards):
    dge_type = wildcards.dge_type

    extra_params = ""

    if dge_type == ".exon":
        extra_params = ""
    elif dge_type == ".intron":
        extra_params = "LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    elif dge_type == ".all":
        extra_params = "LOCUS_FUNCTION_LIST=INTRONIC"
    if dge_type == ".Reads_exon":
        extra_params = "OUTPUT_READS_INSTEAD=true"
    elif dge_type == ".Reads_intron":
        extra_params = (
            "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=null"
            + "LOCUS_FUNCTION_LIST=INTRONIC"
        )
    elif dge_type == ".Reads_all":
        extra_params = "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=INTRONIC"

    if wildcards.mm_included == ".mm_included":
        extra_params = extra_params + " READ_MQ=0"

    return extra_params


def get_files_to_merge(pattern, project_id, sample_id, **kwargs):
    # recursive function to find all files to merge. a merged sample can be merged
    # from merged samples. to avoid cyclic dependencies, here we look for all files
    # which are the dependencies of the underlying samples
    is_merged = project_df.get_metadata(
        "is_merged", project_id=project_id, sample_id=sample_id
    )

    files = []

    if not is_merged:
        files = expand(pattern, sample_id=sample_id, project_id=project_id, **kwargs)
    else:
        merge_ix = project_df.get_metadata(
            "merged_from", sample_id=sample_id, project_id=project_id
        )

        for (p, s) in merge_ix:
            files = files + get_files_to_merge(
                project_id=p, sample_id=s, pattern=pattern, **kwargs
            )

    return list(set(files))


def get_files_to_merge_snakemake(pattern):
    # inner function to be returned
    def get_merged_pattern(wildcards):
        kwargs = {}

        # konvert wildcards to dict
        for key, value in wildcards.items():
            kwargs[key] = value

        files = get_files_to_merge(pattern=pattern, **kwargs)

        return files

    return get_merged_pattern


def get_ribo_depletion_log(wildcards):
    is_merged = project_df.get_metadata(
        "is_merged", sample_id=wildcards.sample_id, project_id=wildcards.project_id
    )

    if is_merged:
        return [merged_ribo_depletion_log]
    else:
        return [ribo_depletion_log]


def get_top_barcodes(wildcards):
    if wildcards.n_beads == "spatial":
        return {"top_barcodes": spatial_barcodes}
    if wildcards.dge_cleaned == "":
        return {"top_barcodes": top_barcodes}
    else:
        return {"top_barcodes": top_barcodes_clean}


def get_parsed_puck_file(wildcards):
    is_spatial = project_df.is_spatial(
        project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )

    if is_spatial:
        return {"puck_file": parsed_spatial_barcodes}
    else:
        return []


def get_dge_from_run_mode(
    project_id, sample_id, run_mode, data_root_type, downsampling_percentage
):
    has_dge = project_df.has_dge(
        project_id=project_id, sample_id=sample_id)

    if not has_dge:
        raise SpacemakeError(
            f'Sample with id (project_id, sample_id)={project_id}, {sample_id})' +
            f' does not have a DGE')

    is_spatial = project_df.is_spatial(project_id=project_id, sample_id=sample_id)

    is_external = project_df.is_external(project_id=project_id, sample_id=sample_id)

    run_mode_variables = project_df.config.get_run_mode(run_mode).variables

    dge_type = ""
    dge_cleaned = ""
    polyA_adapter_trimmed = ""
    mm_included = ""

    # assign wildcards only for internal samples
    if not is_external:
        if run_mode_variables["polyA_adapter_trimming"]:
            polyA_adapter_trimmed = ".polyA_adapter_trimmed"

        if run_mode_variables["count_intronic_reads"]:
            dge_type = ".all"
        else:
            dge_type = ".exon"

        if run_mode_variables["count_mm_reads"]:
            mm_included = ".mm_included"

        if run_mode_variables["clean_dge"]:
            dge_cleaned = ".cleaned"

    if run_mode_variables["mesh_type"] == "hexagon":
        spot_diameter_um = run_mode_variables["mesh_spot_diameter_um"]
        spot_distance_um = "hexagon"
    elif run_mode_variables["mesh_type"] == "circle":
        spot_diameter_um = run_mode_variables["mesh_spot_diameter_um"]
        spot_distance_um = run_mode_variables["mesh_spot_distance_um"]

    external_wildcard = ""

    n_beads = run_mode_variables["n_beads"]

    if is_external:
        n_beads = "external"
        external_wildcard = ".external"

    if is_spatial:
        n_beads = "spatial"

    # select which pattern
    # if sample is not spatial, we simply select the normal, umi_filtered
    # dge, with the top_n barcodes
    # otherwise, if sample is spatial, either we return he whole dge, containing
    # all beads, or the a meshgrid
    if not is_spatial:
        dge_out_pattern = dge_out_h5ad
        dge_out_summary_pattern = dge_out_h5ad_obs
    elif run_mode_variables["mesh_data"]:
        dge_out_pattern = dge_spatial_mesh
        dge_out_summary_pattern = dge_spatial_mesh_obs
    else:
        dge_out_pattern = dge_spatial
        dge_out_summary_pattern = dge_spatial_obs

    out_files_pattern = {"dge_summary": dge_out_summary_pattern, "dge": dge_out_pattern}

    out_files = {
        key: expand(
            pattern,
            project_id=project_id,
            sample_id=sample_id,
            dge_type=dge_type,
            dge_cleaned=dge_cleaned,
            polyA_adapter_trimmed=polyA_adapter_trimmed,
            mm_included=mm_included,
            spot_diameter_um=spot_diameter_um,
            spot_distance_um=spot_distance_um,
            n_beads=n_beads,
            is_external=external_wildcard,
            data_root_type=data_root_type,
            downsampling_percentage=downsampling_percentage,
        )
        for key, pattern in out_files_pattern.items()
    }

    return out_files


def get_qc_sheet_input_files(wildcards):
    # returns star_log, reads_type_out, strand_info
    # first checks the run modes, and returns either polyA_adapter_trimmed, untrimmed
    # or both
    project_id = wildcards.project_id
    sample_id = wildcards.sample_id

    is_merged = project_df.get_metadata(
        "is_merged", project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )

    run_modes = get_run_modes_from_sample(wildcards.project_id, wildcards.sample_id)

    is_polyA_adapter_trimmed = set(
        [x["polyA_adapter_trimming"] for x in run_modes.values()]
    )

    # if sample has both polyA trimmed and untrimmed mapped bam files
    if len(is_polyA_adapter_trimmed) == 2:
        polyA_adapter_trimmed_wildcard = ["", ".polyA_adapter_trimmed"]
    elif True in is_polyA_adapter_trimmed:
        polyA_adapter_trimmed_wildcard = [".polyA_adapter_trimmed"]
    elif False in is_polyA_adapter_trimmed:
        polyA_adapter_trimmed_wildcard = [""]

    extra_args = {
        "sample_id": wildcards.sample_id,
        "project_id": wildcards.project_id,
        "polyA_adapter_trimmed": polyA_adapter_trimmed_wildcard,
    }

    if is_merged:
        star_log_pattern = merged_star_log_file
    else:
        star_log_pattern = star_log_file

    to_return = {
        "star_log": expand(star_log_pattern, **extra_args),
        "reads_type_out": expand(reads_type_out, **extra_args),
        "strand_info": expand(strand_info, **extra_args),
    }

    for run_mode in run_modes:
        run_mode_dge = get_dge_from_run_mode(
            project_id,
            sample_id,
            run_mode,
            data_root_type=wildcards.data_root_type,
            downsampling_percentage=wildcards.downsampling_percentage,
        )

        to_return[f"{run_mode}.dge_summary"] = run_mode_dge["dge_summary"]

    return to_return


def get_bam_tag_names(project_id, sample_id):
    barcode_flavor = project_df.get_metadata(
        "barcode_flavor", project_id=project_id, sample_id=sample_id
    )

    bam_tags = config["barcode_flavors"][barcode_flavor]["bam_tags"]

    tag_names = {}

    for tag in bam_tags.split(","):
        tag_name, tag_variable = tag.split(":")

        tag_names[tag_variable] = tag_name

    return tag_names


def get_puck_file(wildcards):
    if not project_df.is_spatial(
        project_id=wildcards.project_id, sample_id=wildcards.sample_id
    ):

        return []

    puck_barcode_file = project_df.get_metadata(
        "puck_barcode_file", project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )

    puck_name = project_df.get_metadata(
        "puck", project_id=wildcards.project_id, sample_id=wildcards.sample_id
    )
    
    puck = project_df.config.get_puck(puck_name)

    if puck_barcode_file is None:
        return {"barcode_file": puck.variables['barcodes']}
    else:
        return {"barcode_file": puck_barcode_file}


def get_automated_analysis_dge_input(wildcards):
    # there are three options:
    # 1) no spatial dge
    # 2) spatial dge, no mesh
    # 3) spatial dge with a mesh
    return [
        get_dge_from_run_mode(
            project_id=wildcards.project_id,
            sample_id=wildcards.sample_id,
            run_mode=wildcards.run_mode,
            data_root_type=wildcards.data_root_type,
            downsampling_percentage=wildcards.downsampling_percentage,
        )["dge"]
    ]

def get_novosparc_input_files(config):
    if ('reference_project_id' in config and
        config['reference_project_id'] != '' and
        'reference_sample_id' in config and
        config['reference_sample_id'] != '' and
        'reference_run_mode' in config and
        config['reference_run_mode'] != '' and
        'reference_umi_cutoff' in config and
        config['reference_umi_cutoff'] != ''):

        ret = expand(novosparc_with_reference_h5ad,
            data_root_type = 'complete_data',
            downsampling_percentage = '',
            project_id = config['project_id'],
            sample_id = config['sample_id'],
            run_mode = config['run_mode'],
            umi_cutoff = config['umi_cutoff'],
            reference_project_id = config['reference_project_id'],
            reference_sample_id = config['reference_sample_id'],
            reference_run_mode = config['reference_run_mode'],
            reference_umi_cutoff = config['reference_umi_cutoff'],
        )

    elif ('project_id' in config and
          'sample_id' in config and
          'umi_cutoff' in config and
          'run_mode' in config):
        ret = expand(novosparc_denovo_h5ad,
            data_root_type = 'complete_data',
            downsampling_percentage = '',
            project_id = config['project_id'],
            sample_id = config['sample_id'],
            run_mode = config['run_mode'],
            umi_cutoff = config['umi_cutoff'],
        )
    else:
        ret = []
    
    return ret

def get_novosparc_with_reference_input_files(wildcards):
    out_dict = {}

    st_adata = expand(automated_analysis_result_file,
        data_root_type = wildcards.data_root_type,
        downsampling_percentage = wildcards.downsampling_percentage,
        # access data without the first underscore
        project_id = wildcards.reference_project_id,
        sample_id = wildcards.reference_sample_id,
        run_mode = wildcards.reference_run_mode,
        umi_cutoff = wildcards.reference_umi_cutoff,
    )

    out_dict['st_adata'] = st_adata

    return out_dict

