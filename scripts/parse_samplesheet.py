import pandas as pd
from os.path import join

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
    ss = pd.read_csv(fp_samplesheet, sep=",", skiprows=21)

    # bcl2fasta automatically changes - into _ char in output filenames
    for f in ['Sample_ID', 'Sample_Name', 'Sample_Project']:
        ss[f] = ss[f].apply(lambda x: x.replace('-', '_') if type(x) != float else x)

    fp_fastqs = []
    for idx, row in ss.iterrows():
        fp_fastq = ''
        if pd.notnull(row['Sample_Project']):
            fp_fastq = row['Sample_Project']
        if pd.notnull(row['Sample_Name']):
            fp_fastq = join(fp_fastq, row['Sample_ID'])
        for direction in ['R1', 'R2']:
            fp_fastqs.append(join(fp_fastq, '%s_S%i_L%03i_%s_001.fastq.gz' % (row['Sample_Name'], idx, int(row['Lane']), direction)))

    print('\n'.join(fp_fastqs))
    return None
