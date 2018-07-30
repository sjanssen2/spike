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


def get_fastq_filenames(fp_samplesheet):
    """Returns list of fastq.gz filepaths parsed from given sample sheet.

    Parameters
    ----------
    fp_samplesheet : str
        Filepath of sample sheet

    Returns
    -------
    List of filepaths.
    """
    ss = parse_samplesheet(fp_samplesheet)

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


def get_sample_fastqprefixes(fp_samplesheet):
    ss = parse_samplesheet(fp_samplesheet)
    return list(ss['fastq-prefix'].unique())


def get_lanes_for_sampleID(fp_samplesheet, sample):
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
    ss = parse_samplesheet(fp_samplesheet)
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


def get_sample_names(fp_samplesheet, ukd_actions={'trio', 'somatic'}):
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

    ss = parse_samplesheet(fp_samplesheet)
    samples = []
    for n, g in ss[ss['ukd_action'].isin(ukd_actions)].fillna('').groupby(['Sample_Project', 'Sample_Name', 'Sample_ID', 's-idx']):
        samples.append(g['fastq-prefix'].unique()[0])
    return samples


def get_role(ukd_project, ukd_entity_id, ukd_entity_role, config):
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
    # parse all available sample sheets
    fps_samplesheets = glob.glob(join(config['dirs']['prefix'], config['dirs']['inputs'], config['dirs']['samplesheets'], '*XX_ukd.csv'))

    global_samplesheet = []
    for fp_samplesheet in fps_samplesheets:
        ss = parse_samplesheet(fp_samplesheet)
        ss['run'] = fp_samplesheet.split('/')[-1].replace('_ukd.csv', '')
        global_samplesheet.append(ss)
    global_samplesheet = pd.concat(global_samplesheet)

    samples = global_samplesheet

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

    res = {join(sample['run'], sample['fastq-prefix']) for idx, sample in samples.iterrows()}
    if len(res) > 1:
        raise ValueError("Stefan, check if use cases can occour with more than one result!")

    return list(res)[0]
