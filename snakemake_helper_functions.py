def hamming_distance(string1, string2):
    return sum(c1 != c2 for c1, c2 in zip(string1, string2))

def compute_max_barcode_mismatch(indices):
    """computes the maximum number of mismatches allowed for demultiplexing based
    on the indices present in the sample sheet."""
    num_samples = len(indices)
    
    if num_samples == 1:
        return 4
    else:
        max_mismatch = 3
        for i in range(num_samples-1):
            for j in range(i+1, num_samples):
                hd = hamming_distance(indices[i], indices[j])
                max_mismatch = min(max_mismatch, math.ceil(hd/2)-1)
    return max_mismatch

def read_sample_sheet(sample_sheet_path, flowcell_id):
    with open(sample_sheet_path) as sample_sheet:
        ix = 0
        investigator = 'none'
        sequencing_date = 'none'

        for line in sample_sheet:
            line = line.strip('\n')
            if 'Investigator' in line:
                investigator = line.split(',')[1]
            if 'Date' in line:
                sequencing_date = line.split(',')[1]
            if '[Data]' in line:
                break
            else:
                ix = ix + 1

    df = pd.read_csv(sample_sheet_path, skiprows = ix+1)
    df['species'] = df['Description'].str.split('_').str[-1]
    df['investigator'] = investigator
    df['sequencing_date'] = sequencing_date

    # mock R1 and R2
    df['R1'] = 'none'
    df['R2'] = 'none'

    df.rename(columns={"Sample_ID":"sample_id", "Sample_Name":"puck_id", "Sample_Project":"project_id", "Description": "experiment"}, inplace=True)
    
    df['flowcell_id'] = flowcell_id
    df['demux_barcode_mismatch'] = compute_max_barcode_mismatch(df['index'])
    df['sample_sheet'] = sample_sheet_path
    df['demux_dir'] = df['sample_sheet'].str.split('/').str[-1].str.split('.').str[0]

    return df[['sample_id', 'puck_id', 'project_id', 'sample_sheet', 'flowcell_id',
               'species', 'demux_barcode_mismatch', 'demux_dir', 'R1', 'R2', 'investigator', 'sequencing_date', 'experiment']]    

def create_lookup_table(df):
    samples_lookup = {}

    projects = df.project_id.unique()

    for p in projects:
        sample_ids = df[df.project_id.eq(p)].sample_id.to_list()
        sample_sheet = df[df.project_id.eq(p)].sample_sheet.to_list()[0]
        flowcell_id = df[df.project_id.eq(p)].flowcell_id.to_list()[0]
        demux_barcode_mismatch = df[df.project_id.eq(p)].demux_barcode_mismatch.to_list()[0]
        demux_dir = df[df.project_id.eq(p)].demux_dir.to_list()[0]
        investigator = df[df.project_id.eq(p)].investigator.to_list()[0]
        sequencing_date = df[df.project_id.eq(p)].sequencing_date.to_list()[0]
       
        # per sample
        species = df[df.project_id.eq(p)].species.to_list()
        pucks = df[df.project_id.eq(p)].puck_id.to_list()
        R1 = df[df.project_id.eq(p)].R1.to_list()
        R2 = df[df.project_id.eq(p)].R2.to_list()
        experiment = df[df.project_id.eq(p)].experiment.to_list()

        samples = {}
        for i in range(len(sample_ids)):
            samples[sample_ids[i]] = {
                'species': species[i],
                'puck': pucks[i],
                'R1': R1[i],
                'R2': R2[i],
                'experiment': experiment[i]
            }

        samples_lookup[p] = {
            'sample_sheet': sample_sheet,
            'flowcell_id': flowcell_id,
            'samples': samples,
            'demux_barcode_mismatch': int(demux_barcode_mismatch),
            'demux_dir': demux_dir,
            'investigator': investigator,
            'sequencing_date': sequencing_date
        }

    return samples_lookup

def get_demux_dir(wildcards):
    return samples[wildcards.project]['demux_dir']

def get_demux_indicator(wildcards):
    demux_dir = get_demux_dir(wildcards) 

    return expand(demux_indicator, demux_dir=demux_dir)

def get_species_info(wildcards):
    # This function will return 3 things required by STAR:
    #    - annotation (.gtf file)
    #    - genome (.fa file)
    #    - index (a directory where the STAR index is)
    species = samples[wildcards.project]['samples'][wildcards.sample]['species']

    return {
        'annotation': config['knowledge']['annotations'][species],
        'genome': config['knowledge']['genomes'][species],
        'index': config['knowledge']['indices'][species]['star']
    }

def get_dge_extra_params(wildcards):
    dge_type = wildcards.dge_type

    if dge_type == '_exon':
        return ''
    elif dge_type == '_intron':
        return "LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    elif dge_type == '_all':
        return "LOCUS_FUNCTION_LIST=INTRONIC"
    if dge_type == 'Reads_exon':
        return "OUTPUT_READS_INSTEAD=true"
    elif dge_type == 'Reads_intron':
        return "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    elif dge_type == 'Reads_all':
        return "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=INTRONIC"

def get_basecalls_dir(wildcards):
    flowcell_id = samples[demux_dir2project[wildcards.demux_dir]]['flowcell_id']
    
    basecalls_dir = '/data/remote/basecalls/'
    local_nextseq_raw = '/data/rajewsky/sequencing/nextSeqRaw/'

    # check if flowcell_id exists in /data/remote/basecalls
    if os.path.isdir(basecalls_dir + flowcell_id):
        return [basecalls_dir + flowcell_id]
    # if not, check if flowcell_id exists in local /data/rajewsky/sequencing/nextSeqRaw
    elif os.path.isdir(local_nextseq_raw + flowcell_id):
        return [local_nextseq_raw + flowcell_id]
    # else return a fake path, which won't be present, so snakemake will fail for this, as input directory will be missing
    else:
        return [basecalls_dir + 'none'] 

###############################
# Joining optical to illumina #
###############################

def get_sample_info(raw_folder):
    batches = os.listdir(raw_folder)

    df = pd.DataFrame(columns=['batch_id', 'puck_id'])

    for batch in batches:
        batch_dir = microscopy_raw + '/' + batch

        puck_ids = os.listdir(batch_dir)

        df = df.append(pd.DataFrame({'batch_id': batch, 'puck_id': puck_ids}), ignore_index=True)
        
    return df


###################
# Merging samples #
###################
def get_project(sample):
    # return the project id for a given sample id
    return project_df[project_df.sample_id.eq(sample)].project_id.to_list()[0]

def get_dropseq_final_bam(wildcards):
    # merged_name contains all the samples which should be merged,
    # separated by a dot each
    samples = config['samples_to_merge'][wildcards.merged_project][wildcards.merged_sample]

    input_bams = []

    for sample in samples:
        input_bams = input_bams + expand(dropseq_final_bam,
            project = get_project(sample),
            sample = sample)
    return input_bams

def get_merged_bam_inputs(wildcards):
    # currently not used as we do not tag the bam files with the sample name
    samples = config['samples_to_merge'][wildcards.merged_project][wildcards.merged_sample]

    input_bams = []

    for sample in samples:
        input_bams = input_bams + expand(sample_tagged_bam, 
                merged_sample = wildcards.merged_name,
                sample = sample)

    return input_bams

def get_merged_star_log_inputs(wildcards):
    samples = config['samples_to_merge'][wildcards.merged_project][wildcards.merged_sample]
    
    input_logs = []

    for sample in samples:
        input_logs = input_logs + expand(star_log_file,
                project = get_project(sample),
                sample = sample)

    return input_logs
