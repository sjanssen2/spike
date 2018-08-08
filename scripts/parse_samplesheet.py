import yaml
import pandas as pd
import numpy as np
from os.path import join
import glob


DIRECTIONS = ['R1', 'R2']

def parse_samplesheet(fp_samplesheet):
    ss = pd.read_csv(fp_samplesheet, sep=",", skiprows=21, dtype={'Sample_Name': str, 'Sample_ID': str})

    # bcl2fasta automatically changes - into _ char in output filenames
    for f in ['Sample_ID', 'Sample_Name', 'Sample_Project']:
        ss[f] = ss[f].apply(lambda x: x.replace('-', '_') if type(x) != float else x)

    # bcl2fastq uses a S%03i index to address samples.
    # They are numbered as occuring in the samplesheet order starting with 1.
    # However, number is not increased if Sample_ID was already seen.
    uidx = dict()
    for _, sample_id in ss['Sample_ID'].iteritems():
        if sample_id not in uidx:
            uidx[sample_id] = len(uidx) + 1
    ss['s-idx'] = ss['Sample_ID'].apply(lambda x: uidx[x])

    # TODO: ensure that sample names do not clash when not considering s-idx!

    # fastq-prefix
    fp_fastqs = []
    for idx, row in ss.iterrows():
        fp_fastq = ''
        if pd.notnull(row['Sample_Project']):
            fp_fastq = row['Sample_Project']
        if pd.notnull(row['Sample_Name']):
            fp_fastq = join(fp_fastq, row['Sample_ID'])
        fp_fastqs.append(join(fp_fastq,
            '%s' % (
                row['Sample_Name'] if pd.notnull(
                    row['Sample_Name']) else row['Sample_ID'])))
    ss['fastq-prefix'] = fp_fastqs

    return ss


def get_fastq_filenames(run, samplesheets):
    """Returns list of fastq.gz filepaths parsed from given sample sheet.

    Parameters
    ----------
    fp_samplesheet : str
        Filepath of sample sheet

    Returns
    -------
    List of filepaths.
    """
    ss = samplesheets[samplesheets['run'] == run]

    fp_fastqs = []
    for idx, row in ss.iterrows():
        fp_fastq = row['fastq-prefix']
        for direction in DIRECTIONS:
            fp_fastqs.append(
                '%s_L%03i_%s_001.fastq.gz' % (
                    fp_fastq,
                    int(row['Lane']),
                    direction))

    # add fps for undetermined reads
    for lane in ss['Lane'].unique():
        for direction in DIRECTIONS:
            fp_fastqs.append(
                'Undetermined_L%03i_%s_001.fastq.gz' % (lane, direction))

    return fp_fastqs


def get_sample_fastqprefixes(run, samplesheets):
    ss = samplesheets[samplesheets['run'] == run]
    return list(ss['fastq-prefix'].unique())


def get_lanes_for_sampleID(sample, samplesheets):
    """Return lanes a given sample is spread across.

    Parameters
    ----------
    fp_samplesheet : str
        Filepath to Sample Sheet
    sampleName : str
        The name of the sample for which lanes should be looked up.
    sampleID : str
        The ID of the sample for which lanes should be looked up.
    sidx : str
        The running index of the sample for which lanes should be looked up.

    Returns
    -------
    [str] : Lane numbers on which to find the given sample.
    """
    ss = samplesheets.copy()
    ss['tmp-id'] = ['%s/%s%s' % (row['Sample_Project'], row['Sample_ID'], '/'+row['Sample_Name'] if pd.notnull(row['Sample_Name']) else "") for _, row in ss.iterrows()]
    res = ss[ss['tmp-id'] == sample]['Lane'].unique()

    return res


def get_laneSplitInputs(wildcards, dir_input_samplesheets, dir_intermediate_demultiplex):
    """Given targeted joined output fastq.gz file, obtain fastq.gz files eventually split across lanes."""
    #params must contain at least
    #{'prefix': '/home/jansses/gpfs/', 'run': '180614_SN737_0438_BCC7MCACXX', 'sample': 'Alps/ALPS_66_a_S18', 'direction': 'R1'}
    params = dict(wildcards)
    fp_samplesheet = join(params['prefix'], dir_input_samplesheets, params['run']) + '_ukd.csv'
    print("params", params)
    print("fp_samplesheet", fp_samplesheet)
    ss = parse_samplesheet(fp_samplesheet)
    ss['tmp-id'] = ['%s/%s%s_S%s' % (row['Sample_Project'], row['Sample_ID'], '/'+row['Sample_Name'] if pd.notnull(row['Sample_Name']) else "", row['s-idx']) for _, row in ss.iterrows()]

    lanes = ss[ss['tmp-id'] == params['sample']]['Lane'].unique()

    res = ["%s_L%03i_%s_001.fastq.gz" % (
        join(params['prefix'], dir_intermediate_demultiplex, params['run'], params['sample']),
        int(lane),
        params['direction']) for lane in lanes]
    # print(params, res, ss['tmp-id'])
    return res


def get_sample_names(samplesheets, ukd_actions={'trio', 'somatic'}):
    """Get unique list of full sample names, i.e. with Project-Name, Sample-Group-Name and s-idx.

    Parameters
    ----------
    fp_samplesheet : str
        Filepath to SampleSheet
    ukd_action : set(str)
        Default: {'trio', 'somatic'}.
        Limit returned samples to those having one of the ukd_actions assigned

    Returns
    -------
    [str] the unique full qualified sample names."""

    ss = samplesheets
    samples = []
    for n, g in ss[ss['ukd_action'].isin(ukd_actions)].fillna('').groupby(['Sample_Project', 'Sample_Name', 'Sample_ID', 's-idx']):
        samples.append(g['fastq-prefix'].unique()[0])
    return samples


def get_global_samplesheets(dir_samplesheets):
    # parse all available sample sheets
    fps_samplesheets = glob.glob('%s*XX_ukd.csv' % dir_samplesheets)

    global_samplesheet = []
    for fp_samplesheet in fps_samplesheets:
        ss = parse_samplesheet(fp_samplesheet)
        ss['run'] = fp_samplesheet.split('/')[-1].replace('_ukd.csv', '')
        global_samplesheet.append(ss)
    global_samplesheet = pd.concat(global_samplesheet, sort=False)

    return global_samplesheet


def get_role(ukd_project, ukd_entity_id, ukd_entity_role, samplesheets):
    """Returns file path for bam, given project, entity and role (for trio).

    Parameters
    ----------
    ukd_project : str
        Name of project, to avoid entity ID clashes across projects.
    ukd_entity_id : str
        Entity ID for which role needs to be obtained.
    ukd_entity_role : str
        Role of entity ID whose bam filepath shall be returned.
    config : snakemake.config
        Snakemakes config object to obtain file path of sample sheets.

    Returns
    -------
    str: Filepath of bam file for given entity role.
    """
    samples = samplesheets

    # select correct project
    try:
        x = samples[samples['Sample_Project'] == ukd_project]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find an UKD project with name "%s". Available projects are:\n\t%s\n' % (ukd_project, '\n\t'.join(sorted(samples['Sample_Project'].unique()))))
    else:
        samples = x

    # select correct entity
    try:
        x = samples[samples['ukd_entity_id'] == ukd_entity_id]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find an UKD entity group with name "%s". Available entities for projects "%s" are:\n\t%s\n' % (ukd_entity_id, ukd_project, '\n\t'.join(sorted(samples['ukd_entity_id'].unique()))))
    else:
        samples = x

    # select correct role
    try:
        x = samples[samples['ukd_entity_role'] == ukd_entity_role]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find a role "%s" for UKD entity group with name "%s". Available roles are:\n\t%s\n' % (ukd_entity_role, ukd_entity_id, '\n\t'.join(sorted(samples['ukd_entity_role'].unique()))))
    else:
        samples = x

    res = {sample['fastq-prefix'] for idx, sample in samples.iterrows()}
    if len(res) > 1:
        raise ValueError("Stefan, check if use cases can occour with more than one result!")

    return list(res)[0]


def get_xenograft_host(fp_samplesheet, sample, samplesheets, config):
    """Returns xenograft host genome name for given sample.

    Parameters
    ----------
    fp_samplesheet : str
        Filepath to sample sheet.
    sample : str
        Fastq prefix for sample, e.g. Alps/ALPS_66_a

    Returns
    -------
    str : empty string if no xenograft, otherwise name of host species name.

    Raises
    ------
    ValueError if:
        a) sample sheet lists conflicting host species names.
        b) an unknown xenograft host species is listed in sample sheet.
    """
    ss = samplesheets

    host_species = list(ss[ss['fastq-prefix'] == sample]['ukd_xenograft_species'].dropna().unique())
    if len(host_species) > 1:
        raise ValueError("Ambiguous xenograft host species in sample sheet!")
    if len(host_species) == 0:
        return ""

    res = ""
    try:
        res += config['xenograft']['references']['homo sapiens']
        res += '_' + config['xenograft']['references'][host_species[0]]
    except KeyError:
        raise ValueError("Unknown xenograft host species!")

    return res


def get_species(sample, samplesheets, config):
    # sample can be a single sample ...
    projects = samplesheets[samplesheets['fastq-prefix'] == sample]['Sample_Project'].unique()

    # ... or an entity
    if len(projects) == 0:
        projects = samplesheets[(samplesheets['Sample_Project'] == sample.split('/')[0]) & (samplesheets['ukd_entity_id'] == sample.split('/')[-1])]['Sample_Project'].unique()

    if len(projects) > 1:
        raise ValueError("Ambiguous projects: '%s' for sample '%s'" % (projects, sample))

    return config['projects'][projects[0]]['species']


def get_reference_genome(sample, samplesheets, config):
    return config['references']['genomes'][get_species(sample, samplesheets, config)]

def get_reference_knowns(sample, samplesheets, config, _key):
    return [k for k in config['references']['knowns'][get_species(sample, samplesheets, config)] if _key in k]

def get_reference_exometrack(sample, samplesheets, config):
    return config['references']['exometrack'][get_species(sample, samplesheets, config)]['file']

def get_reference_varscan_somatic(sample, samplesheets, config):
    return config['references']['varscan_somatic'][get_species(sample, samplesheets, config)]

def get_reference_DBSNP(sample, samplesheets, config):
    res = config['references']['DBSNP'][get_species(sample, samplesheets, config)]
    if res is None:
        res = []
    return res


######## avoid run
def _run2date(run):
    date='%04i/%02i/%02i' % (
        int(run.split('_')[0][:2])+2000,
        int(run.split('_')[0][3:4]),
        int(run.split('_')[0][5:6]))
    return date

def get_bwa_mem_header(sample, samplesheets, config):
    samples = samplesheets[samplesheets['fastq-prefix'] == sample]
    res = ' -R "@RG\\tID:%s\\tCN:Department_of_Pediatric_Oncology_Dusseldorf\\tPU:%s\\tDT:%s\\tPL:ILLUMINA\\tLB:%s\\tSM:readgroups.info"' % (
        ' and '.join(samples['run'].dropna().unique()),
        ' and '.join(list(map(lambda x: x.split('_')[-1][1:], samples['run'].dropna().unique()))),
        ' and '.join(list(map(_run2date, samples['run'].dropna().unique()))),
        config['references']['exometrack'][get_species(sample, samplesheets, config)]['protocol_name']
        )
    return res


def get_demux_samples(samplesheets, config):
    # get projects that require snv vs. reference analysis
    background_projects = [prj_name for prj_name in config['projects'] if 'demultiplex' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to tumor vs. normal projects
    background_samples = samplesheets[samplesheets['Sample_Project'].isin(background_projects)]

    return list(background_samples['run'].unique())
    # samples = []
    # for _, sample in background_samples.iterrows():
    #     samples.extend(['%s/%s_L%03i_%s_001.fastq.gz' % (sample['run'], sample['fastq-prefix'], sample['Lane'], direction) for direction in config['directions']])
    #
    # return samples


def get_samples(samplesheets, config):
    # get projects that require snv vs. reference analysis
    background_projects = [prj_name for prj_name in config['projects'] if 'background' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to tumor vs. normal projects
    background_samples = samplesheets[samplesheets['Sample_Project'].isin(background_projects)]

    samples = []
    for sample, g in background_samples.groupby(['Sample_Project', 'fastq-prefix']):
        samples.append({'Sample_Project': sample[0],
                        'sample': sample[1]})

    return samples


def get_tumorNormalPairs(samplesheets, config):
    # get projects that require tumor vs. normal analysis
    tumornormal_projects = [prj_name for prj_name in config['projects'] if 'tumornormal' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to tumor vs. normal projects
    tumornormal_samples = samplesheets[samplesheets['Sample_Project'].isin(tumornormal_projects)]

    pairs = []
    for pair, g in tumornormal_samples.groupby(['Sample_Project', 'ukd_entity_id']):
        # only choose comlete pairs
        if set(g['ukd_entity_role'].unique()) == {'healthy','tumor'}:
            pairs.append({'Sample_Project': pair[0],
                          'ukd_entity_id': pair[1]})

    return pairs


def get_trios(samplesheets, config):
    # get projects that require trio analysis
    trio_projects = [prj_name for prj_name in config['projects'] if 'trio' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to trio projects
    trio_samples = samplesheets[samplesheets['Sample_Project'].isin(trio_projects)]

    trios = []
    for trio, g in trio_samples.groupby(['Sample_Project', 'ukd_entity_id']):
        # only choose comlete trios
        if set(g['ukd_entity_role'].unique()) == {'patient', 'mother', 'father'}:
            trios.append({'Sample_Project': trio[0],
                          'ukd_entity_id': trio[1]})

    return trios


def get_projects_with_exomecoverage(config):
    res = []
    for name in config['projects']:
        if ('actions' in config['projects'][name]) and (len(set(config['projects'][name]['actions']) & set(['background', 'trio', 'tumornormal'])) > 0):
            res.append(name)
    return res


def get_rejoin_fastqs(sample, samplesheets, config):
    res = []
    for _, row in samplesheets[samplesheets['fastq-prefix'] == sample].iterrows():
        res.append('%s/%s_L%03i' % (row['run'], row['fastq-prefix'], row['Lane']))
    return res
