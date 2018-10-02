import yaml
import pandas as pd
import numpy as np
from os.path import join
import glob
import sys


def parse_samplesheet(fp_samplesheet):
    ss = pd.read_csv(fp_samplesheet, sep=",", skiprows=21, dtype={'Sample_Name': str, 'Sample_ID': str})

    # bcl2fasta automatically changes - into _ char in output filenames
    idx_rawilluminainput = ss[pd.notnull(ss['Lane'])].index
    for f in ['Sample_ID', 'Sample_Name', 'Sample_Project']:
        ss.loc[idx_rawilluminainput, f] = ss.loc[idx_rawilluminainput, f].apply(lambda x: x.replace('-', '_') if type(x) != float else x)

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

    # remove samples that are marked to be ignored
    ss = ss[pd.isnull(ss['spike_ignore_sample'])]

    # set Lane to 0 if not defined, as in Macrogen samples
    ss['Lane'] = ss['Lane'].fillna(0)

    return ss


def get_global_samplesheets(dir_samplesheets):
    # parse all available sample sheets
    fps_samplesheets = glob.glob('%s*_spike.csv' % dir_samplesheets)

    global_samplesheet = []
    for fp_samplesheet in fps_samplesheets:
        ss = parse_samplesheet(fp_samplesheet)
        ss['run'] = fp_samplesheet.split('/')[-1].replace('_spike.csv', '')
        global_samplesheet.append(ss)
    global_samplesheet = pd.concat(global_samplesheet, sort=False)

    return global_samplesheet


def get_role(spike_project, spike_entity_id, spike_entity_role, samplesheets):
    """Returns file path for bam, given project, entity and role (for trio).

    Parameters
    ----------
    spike_project : str
        Name of project, to avoid entity ID clashes across projects.
    spike_entity_id : str
        Entity ID for which role needs to be obtained.
    spike_entity_role : str
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
        x = samples[samples['Sample_Project'] == spike_project]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find an spike project with name "%s". Available projects are:\n\t%s\n' % (spike_project, '\n\t'.join(sorted(samples['Sample_Project'].unique()))))
    else:
        samples = x

    # select correct entity
    try:
        x = samples[samples['spike_entity_id'] == spike_entity_id]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find an spike entity group with name "%s". Available entities for projects "%s" are:\n\t%s\n' % (spike_entity_id, spike_project, '\n\t'.join(sorted(samples['spike_entity_id'].unique()))))
    else:
        samples = x

    # select correct role
    try:
        x = samples[samples['spike_entity_role'] == spike_entity_role]
        x.iloc[0]
    except IndexError:
        raise ValueError('Could not find a role "%s" for spike entity group with name "%s". Available roles are:\n\t%s\n' % (spike_entity_role, spike_entity_id, '\n\t'.join(sorted(samples['spike_entity_role'].unique()))))
    else:
        samples = x

    res = {sample['fastq-prefix'] for idx, sample in samples.iterrows()}
    if len(res) > 1:
        raise ValueError("Stefan, check if use cases can occour with more than one result!")

    return list(res)[0]


def get_species(sample, samplesheets, config):
    # sample can be a single sample ...
    projects = samplesheets[samplesheets['fastq-prefix'] == sample]['Sample_Project'].unique()

    # ... or an entity
    if len(projects) == 0:
        projects = samplesheets[(samplesheets['Sample_Project'] == sample.split('/')[0]) & (samplesheets['spike_entity_id'] == sample.split('/')[-1])]['Sample_Project'].unique()

    if len(projects) > 1:
        raise ValueError("Ambiguous projects: '%s' for sample '%s'" % (projects, sample))

    return config['projects'][projects[0]]['species']


def get_reference_genome(sample, samplesheets, config):
    return config['references']['genomes'][get_species(sample, samplesheets, config)]


def get_reference_knowns(sample, samplesheets, config, _key):
    return [k for k in config['references']['knowns'][get_species(sample, samplesheets, config)] if _key in k]


def get_reference_exometrack(sample, samplesheets, config, returnfield='file'):
    # there might be a project specific exome track, like for samples we got sequenced by macrogen:
    projects = samplesheets[
        (samplesheets['fastq-prefix'] == sample) | # check samplenames
        (('%s/%s' % (samplesheets['Sample_Project'], samplesheets['spike_entity_id'])) == sample) # or project-name / spike_entity
        ]['Sample_Project'].unique()
    if len(projects) != 1:
        raise ValueError('Ambigious or missing project for sample "%s"!' % sample)
    if 'exometrack' in config['projects'][projects[0]]:
        return config['projects'][projects[0]]['exometrack'][returnfield]

    # by default, we return the species specific exome track
    return config['references']['exometrack'][get_species(sample, samplesheets, config)][returnfield]


def get_reference_varscan_somatic(sample, samplesheets, config):
    return config['references']['varscan_somatic'][get_species(sample, samplesheets, config)]


######## avoid run
def _run2date(run):
    date='%04i/%02i/%02i' % (
        int(run.split('_')[0][:2])+2000,
        int(run.split('_')[0][3:4]),
        int(run.split('_')[0][5:6]))
    return date


def get_bwa_mem_header(sample, samplesheets, config):
    samples = samplesheets[samplesheets['fastq-prefix'] == sample]
    res = ' -R "@RG\\tID:%s\\tCN:Department_of_Pediatric_Oncology_Dusseldorf\\tPU:%s\\tDT:%s\\tPL:ILLUMINA\\tLB:%s\\tSM:%s"' % (
        ' and '.join(samples['run'].dropna().unique()),
        ' and '.join(list(map(lambda x: x.split('_')[-1][1:], samples['run'].dropna().unique()))),
        ' and '.join(list(map(_run2date, samples['run'].dropna().unique()))),
        get_reference_exometrack(sample, samplesheets, config, returnfield='protocol_name'),
        ' and '.join(samples['Sample_ID'].dropna().unique()),
        )
    return res


def get_demux_samples(samplesheets, config):
    # get projects that require snv vs. reference analysis
    background_projects = [prj_name for prj_name in config['projects'] if 'demultiplex' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to tumor vs. normal projects
    background_samples = samplesheets[samplesheets['Sample_Project'].isin(background_projects)]

    return list(background_samples['run'].unique())


def get_persamplefastq_samples(samplesheets, config, get='runs'):
    # get projects that require snv vs. reference analysis
    background_projects = [prj_name for prj_name in config['projects'] if 'persamplefastq' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to per sample fastq projects
    background_samples = samplesheets[samplesheets['Sample_Project'].isin(background_projects)]

    if get=='samples':
        return list(background_samples['fastq-prefix'].unique())
    elif get=='projects':
        return background_projects

    return list(background_samples['run'].unique())


def get_samples(samplesheets, config):
    # get projects that require snv vs. reference analysis
    background_projects = [prj_name for prj_name in config['projects'] if 'background' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to tumor vs. normal projects
    background_samples = samplesheets[samplesheets['Sample_Project'].isin(background_projects)]

    samples = []
    for sample, g in background_samples.groupby(['Sample_Project', 'fastq-prefix']):
        samples.append({'Sample_Project': sample[0],
                        'sample': sample[1],
                        'spike_entity_id': g['spike_entity_id'].iloc[0]})

    return samples


def get_tumorNormalPairs(samplesheets, config, species=None):
    # get projects that require tumor vs. normal analysis
    tumornormal_projects = [prj_name for prj_name in config['projects'] if 'tumornormal' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to tumor vs. normal projects
    tumornormal_samples = samplesheets[samplesheets['Sample_Project'].isin(tumornormal_projects)]

    pairs = []
    for pair, g in tumornormal_samples.groupby(['Sample_Project', 'spike_entity_id']):
        # only choose comlete pairs
        if set(g['spike_entity_role'].unique()) == {'healthy','tumor'}:
            if species is not None:
                if get_species(g['fastq-prefix'].iloc[0], samplesheets, config) != species:
                    continue
            pairs.append({'Sample_Project': pair[0],
                          'spike_entity_id': pair[1]})

    return pairs


def get_trios(samplesheets, config):
    # get projects that require trio analysis
    trio_projects = [prj_name for prj_name in config['projects'] if 'trio' in config['projects'][prj_name]['actions']]

    # filter samples to those belonging to trio projects
    trio_samples = samplesheets[samplesheets['Sample_Project'].isin(trio_projects)]

    trios = []
    for trio, g in trio_samples.groupby(['Sample_Project', 'spike_entity_id']):
        # only choose comlete trios, i.e. entities that have at least patient, mother and father
        # there might also be siblings, but we ignore those samples for now.
        if len(set(g['spike_entity_role'].unique()) & {'patient', 'mother', 'father'}) == 3:
            trios.append({'Sample_Project': trio[0],
                          'spike_entity_id': trio[1]})

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


def get_xenograft_hybridreference(sample, samplesheets, config):
    project = samplesheets[samplesheets['fastq-prefix'] == sample]['Sample_Project'].dropna().unique()
    if len(project) != 1:
        raise ValueError("_get_reference: Sample '%s' has ambiguous or missing project!" % sample)
    return config['projects'][project[0]]['xenograft']


def get_xenograft_stepname(sample, samplesheets, config):
    project = samplesheets[samplesheets['fastq-prefix'] == sample]['Sample_Project'].dropna().unique()
    if len(project) != 1:
        sys.stderr.write('%s\n' % (project))
        raise ValueError("_get_stepname: Sample '%s' has ambiguous or missing project!" % sample)
    if 'xenograft' in config['projects'][project[0]] and config['projects'][project[0]]['xenograft'] != "":
        return config['stepnames']['xenograft_bwa_sampe']
    else:
        return config['stepnames']['trim']
