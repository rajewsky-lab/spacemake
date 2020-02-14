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
    
    return ['/data/remote/basecalls/' + flowcell_id]
