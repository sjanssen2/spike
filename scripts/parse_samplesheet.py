import pandas as pd
import numpy as np
from os.path import join

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

    # fastq-prefix
    fp_fastqs = []
    for idx, row in ss.iterrows():
        fp_fastq = ''
        if pd.notnull(row['Sample_Project']):
            fp_fastq = row['Sample_Project']
        if pd.notnull(row['Sample_Name']):
            fp_fastq = join(fp_fastq, row['Sample_ID'])
        fp_fastqs.append(join(fp_fastq,
            '%s_S%i' % (
                row['Sample_Name'] if pd.notnull(
                    row['Sample_Name']) else row['Sample_ID'],
                row['s-idx'])))
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
    DIRECTIONS = ['R1', 'R2']

    ss = parse_samplesheet(fp_samplesheet)

    fp_fastqs = []
    for idx, row in ss.iterrows():
        fp_fastq = row['fastq-prefix']
        for direction in DIRECTIONS:
            fp_fastqs.append(
                '_L%s%03i_%s_001.fastq.gz' % (
                    fp_fastq,
                    int(row['Lane']),
                    direction))

    # add fps for undetermined reads
    for lane in ss['Lane'].unique():
        for direction in DIRECTIONS:
            fp_fastqs.append(
                'Undetermined_S0_L%03i_%s_001.fastq.gz' % (lane, direction))

    return fp_fastqs


# def get_sample_fastqprefixes(fp_samplesheet):
#     ss = parse_samplesheet(fp_samplesheet)
#
#     return list(ss['fastq-prefix'].unique())
def get_lanes_for_sampleID(fp_samplesheet, sampleName, sampleID, sidx):
    ss = parse_samplesheet(fp_samplesheet)

    ss['tmp-id'] = ['%s%s%s' % (row['Sample_ID'], '/'+row['Sample_Name'] if pd.notnull(row['Sample_Name']) else "", row['s-idx']) for _, row in ss.iterrows()]
    res = ss[ss['tmp-id'] == '%s%s%s' % (sampleName, sampleID, sidx)]['Lane'].unique()

    return res
