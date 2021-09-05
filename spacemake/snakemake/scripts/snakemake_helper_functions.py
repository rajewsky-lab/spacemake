from spacemake.errors import BarcodeFlavorNotFoundError, SpeciesNotFoundError

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
    for flavor, flavor_settings in config["knowledge"]['barcode_flavor'].items():
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
    flavor = project_df.get_metadata('barcode_flavor', project_id = wildcards.project,
            sample_id = wildcards.sample)
    if flavor not in bc_flavor_data.preprocess_settings:
        raise BarcodeFlavorNotFoundError(flavor)

    settings = bc_flavor_data.preprocess_settings[flavor]

    return settings

def get_demux_indicator(wildcards):
    demux_dir = project_df.get_metadata(
        "demux_dir", sample_id=wildcards.sample, project_id=wildcards.project
    )

    return expand(demux_indicator, demux_dir=demux_dir)


def get_star_input_bam(wildcards):
    if wildcards.polyA_adapter_trimmed == '.polyA_adapter_trimmed':
        return {'reads': tagged_polyA_adapter_trimmed_bam}
    else:
        return {'reads': tagged_or_merged_bam(wildcards)}

def get_mapped_final_bam(wildcards):
    if wildcards.mm_included == '.mm_included':
        return {'reads': final_bam_mm_included_pipe}
    else:
        return {'reads': final_bam}


def get_species_genome_annotation(wildcards):
    # This function will return 2 things required by STAR:
    #    - annotation (.gtf file)
    #    - genome (.fa file)
    if 'species' not in wildcards.keys():
        species = project_df.get_metadata(
            "species", project_id=wildcards.project, sample_id=wildcards.sample
        )
    else:
        species = wildcards.species

    if 'annotations' not in config['knowledge'].keys() or \
            'genomes' not in config['knowledge'].keys() or \
            species not in config['knowledge']['annotations'].keys() or \
            species not in config['knowledge']['genomes'].keys():
        raise SpeciesNotFoundError(species)

    files = {
        "annotation": config["knowledge"]["annotations"][species],
        "genome": config["knowledge"]["genomes"][species]
    }

    return files

def get_star_index(wildcards):
    # This function will return 1 things required by STAR:
    #    - index directory
    species = project_df.get_metadata(
        "species", project_id=wildcards.project, sample_id=wildcards.sample
    )
    return {'index': expand(star_index, species = species)[0]}

def get_rRNA_genome(wildcards):
    return [config['knowledge']['rRNA_genomes'][wildcards.species]]

def get_bt2_rRNA_index(wildcards):
    species = project_df.get_metadata(
        "species", project_id=wildcards.project, sample_id=wildcards.sample
    )

    if 'rRNA_genomes' in config['knowledge']:
        if species in config['knowledge']['rRNA_genomes']:
            return {'index': expand(bt2_rRNA_index_dir, species = species)[0]}
    
    return []

def get_run_mode_variables(run_mode):
    # return the run mode variables
    # first set the default, for each
    # then update each if there is no default
    # load the default
    run_mode_variables = dict(config['run_modes']['default'])

    # first update the default with parent
    if 'parent_run_mode' in config['run_modes'][run_mode].keys():
        parent_run_mode = config['run_modes'][run_mode]['parent_run_mode']
        run_mode_variables.update(config['run_modes'][parent_run_mode])

    run_mode_variables.update(config['run_modes'][run_mode])

    # set the right types
    run_mode_variables['n_beads'] = int(run_mode_variables['n_beads'])
    run_mode_variables['umi_cutoff'] = [int(x) for x  in run_mode_variables['umi_cutoff']]
    run_mode_variables['clean_dge'] = bool(run_mode_variables['clean_dge'])
    run_mode_variables['plot_bead_size'] = float(run_mode_variables['plot_bead_size'])
    run_mode_variables['detect_tissue'] = bool(run_mode_variables['detect_tissue'])
    run_mode_variables['polyA_adapter_trimming'] = bool(run_mode_variables['polyA_adapter_trimming'])
    run_mode_variables['count_mm_reads'] = bool(run_mode_variables['count_mm_reads'])
    run_mode_variables['count_intronic_reads'] = bool(run_mode_variables['count_intronic_reads'])

    return run_mode_variables

def get_run_modes_from_sample(project_id, sample_id):
    run_mode_names = project_df.get_metadata('run_mode', project_id=project_id, sample_id=sample_id)
    
    run_modes = {}

    for run_mode in run_mode_names:
        run_modes[run_mode] = get_run_mode_variables(run_mode)

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
        extra_params = "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=null"+\
                "LOCUS_FUNCTION_LIST=INTRONIC"
    elif dge_type == ".Reads_all":
        extra_params = "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=INTRONIC"

    if wildcards.mm_included == '.mm_included':
        extra_params = extra_params + " READ_MQ=0"

    return extra_params

def get_merged_raw_r2_input(wildcards):
    merge_ix = project_df.get_metadata('merged_from',
        sample_id = wildcards.sample,
        project_id = wildcards.project)

    return [
        expand(raw_reads_mate_2, project = p, sample = s)[0] \ 
        for (p, s) in merge_ix
    ]

def get_merged_bam_input(wildcards):
    merge_ix = project_df.get_metadata('merged_from',
        sample_id = wildcards.sample,
        project_id = wildcards.project)

    return [
        expand(tagged_bam, project = p, sample = s)[0] \
        for (p, s) in merge_ix
    ]

def tagged_or_merged_bam(wildcards):
    is_merged = project_df.get_metadata('is_merged',
        sample_id = wildcards.sample,
        project_id = wildcards.project)

    if is_merged:
        return {'tagged_bam': merged_bam}
    else:
        return {'tagged_bam': tagged_bam}

def get_rRNA_reads_input(wildcards):
    is_merged = project_df.get_metadata('is_merged',
        sample_id = wildcards.sample,
        project_id = wildcards.project)

    if is_merged:
        return {'reads': merged_raw_reads_mate_2}
    else:
        return {'reads': raw_reads_mate_2}

def get_bt2_index(wildcards):
    species = project_df.get_metadata(
        "species", project_id=wildcards.project, sample_id=wildcards.sample
    )

    return config["knowledge"]["indices"][species]["bt2"]

def get_top_barcodes(wildcards):
    if wildcards.dge_cleaned == "":
        return {"top_barcodes": top_barcodes}
    else:
        return {'top_barcodes': top_barcodes_clean}

def get_dge_from_run_mode(
        project_id,
        sample_id,
        run_mode,
        dge_out_pattern = dge_out,
        dge_out_summary_pattern = dge_out_summary,
        **kwargs
    ):
    run_mode_variables = get_run_mode_variables(run_mode)
    
    dge_type = '.exon'
    dge_cleaned = ''
    polyA_adapter_trimmed = ''
    mm_included = ''

    if run_mode_variables['polyA_adapter_trimming']:
        polyA_adapter_trimmed = '.polyA_adapter_trimmed'

    if run_mode_variables['count_intronic_reads']:
        dge_type = '.all'

    if run_mode_variables['count_mm_reads']:
        mm_included = '.mm_included'

    if run_mode_variables['clean_dge']:
        dge_cleaned = '.cleaned'


    dge_out_file = expand(dge_out_pattern,
            project = project_id,
            sample = sample_id,
            dge_type = dge_type,
            dge_cleaned = dge_cleaned,
            polyA_adapter_trimmed = polyA_adapter_trimmed,
            mm_included = mm_included,
            n_beads = run_mode_variables['n_beads'],
            **kwargs)[0]

    dge_out_summary_file = expand(dge_out_summary_pattern,
            project = project_id,
            sample = sample_id,
            dge_type = dge_type,
            dge_cleaned = dge_cleaned,
            polyA_adapter_trimmed = polyA_adapter_trimmed,
            mm_included = mm_included,
            n_beads = run_mode_variables['n_beads'],
            **kwargs)[0]

    return {'dge_summary': dge_out_summary_file,
            'dge': dge_out_file}


def get_dges_from_project_sample(
        project_id,
        sample_id,
        dge_out_pattern = dge_out,
        dge_out_summary_pattern = dge_out_summary,
        **kwargs
    ):
    run_modes = get_run_modes_from_sample(project_id, sample_id).keys()
    dges = {}

    for run_mode in run_modes:
        run_mode_dge = get_dge_from_run_mode(project_id, sample_id, run_mode,
            dge_out_pattern, dge_out_summary_pattern, **kwargs)

        dges[f'{run_mode}.dge'] = run_mode_dge['dge']
        dges[f'{run_mode}.dge_summary'] = run_mode_dge['dge_summary']

    return dges


def get_dge_type(wildcards):
    # expects wildcards to have either a run_mode set, in which case
    # returns one dge+summary pair
    # or project_id and sample_id set and then return all dges
    # associated to this sample (to the run modes of it)
    wildcards_keys = wildcards.keys()
    if 'run_mode' in wildcards_keys:
        return get_dge_from_run_mode(wildcards.project, wildcards.sample,
                wildcards.run_mode)
    elif 'project' in wildcards_keys and 'sample' in wildcards_keys:
        return get_dges_from_project_sample(project_id = wildcards.project,
            sample_id = wildcards.sample)

def get_qc_sheet_input_files(wildcards):
    # returns star_log, reads_type_out, strand_info
    # first checks the run modes, and returns either polyA_adapter_trimmed, untrimmed
    # or both
    run_modes = get_run_modes_from_sample(wildcards.project, wildcards.sample)

    is_polyA_adapter_trimmed = set([x['polyA_adapter_trimming'] for x in run_modes.values()])

    # if sample has both polyA trimmed and untrimmed mapped bam files
    if len(is_polyA_adapter_trimmed) == 2:
        polyA_adapter_trimmed_wildcard = ['', '.polyA_adapter_trimmed']
    elif True in is_polyA_adapter_trimmed:
        polyA_adapter_trimmed_wildcard = ['.polyA_adapter_trimmed']
    elif False in is_polyA_adapter_trimmed:
        polyA_adapter_trimmed_wildcard = ['']

    extra_args = {'sample': wildcards.sample,
                  'project': wildcards.project,
                  'polyA_adapter_trimmed': polyA_adapter_trimmed_wildcard}

    return {
        'star_log': expand(star_log_file, **extra_args),
        'reads_type_out': expand(reads_type_out, **extra_args),
        'strand_info': expand(strand_info, **extra_args)}

def get_bam_tag_names(project_id, sample_id):
    barcode_flavor = project_df.get_metadata('barcode_flavor', project_id = project_id,
            sample_id = sample_id)

    bam_tags = config["knowledge"]["barcode_flavor"][barcode_flavor]["bam_tags"]

    tag_names = {}

    for tag in bam_tags.split(","):
        tag_name, tag_variable = tag.split(":")

        tag_names[tag_variable] = tag_name

    return tag_names

def get_puck_file(wildcards):
    puck_barcode_file = project_df.get_metadata('puck_barcode_file',
            project_id = wildcards.project,
            sample_id = wildcards.sample)

    if puck_barcode_file == "none":
        return []
    else:
        return {"barcode_file": puck_barcode_file}
