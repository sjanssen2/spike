import biom
from biom.util import biom_open
import pandas as pd
import sys


def exclude_sample(sample, blacklist):
    for pattern in blacklist:
        if pattern in sample:
            return True
    return False


def sample_to_biom(fp_gatk_sample, fp_biom, sample=None):
    """Loads a GATK vcf file, extracts SNPs as features and stores as biom file.

    Parameters
    ----------
    fp_gatk_sample : str
        Input file:
        Filepath to 120_gatk_CombineVariants resulting VCF file.
    fp_biom : str
        Output file, resulting BIOM table.
    sample : str
        Default: None
        Spike sample name, e.g. "Keimbahn/KB0056_c"
    """
    with open(fp_gatk_sample, 'r') as f:
        snps = dict()
        for i, line in enumerate(f.readlines()):
            if line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            chromosome, startposition, depths = parts[0], parts[1], None
            for header, value in zip(parts[-2].split(':'), parts[-1].split(':')):
                # Allelic depths for the ref and alt alleles in the order listed
                if header == 'AD':
                    depths = value
                    break
            # currently, feature include only chromosome and start position. One might think about adding ref / alt allel
            feature_id = '-'.join([chromosome, startposition])
            # currently, features are strictly presence / absence, one could use the allelic depth
            snps[feature_id] = 1

    # save feature vector as biom file
    with biom_open(fp_biom, 'w') as f:
        # first, convert data into pandas DataFrame
        if sample is None:
            sample = fp_gatk_sample.split('/')[-1].split('.')[0]
        t = pd.DataFrame(data=None, index=snps.keys(), columns=[sample]).fillna(1)

        # then, store as biom file
        biom.Table(data=t.values,
                   observation_ids=t.index,
                   sample_ids=t.columns
                  ).to_hdf5(f, sample)


def merge_samples(fp_inputs, fp_output, project_name):
    """Take a list of filepaths to biom tables, merge and write to new file.

    Parameters
    ----------
    fp_inputs : [str]
        List of file paths to biom tables.
    fp_output : str
        Filepath of resulting biom table that concats all input biom tables.
    project_name : str
        Name for resulting, merged biom table.
    """
    samples = []
    sys.stderr.write("loading %i biom tables: " % len(fp_inputs))
    for fp_input in fp_inputs:
        sys.stderr.write('.')
        samples.append(biom.load_table(fp_input))
    sys.stderr.write(' done.\nMerging biom tables ...')
    with biom_open(fp_output, 'w') as f:
        samples[0].concat(samples[1:]).to_hdf5(f, project_name)
    sys.stderr.write(' done.\n')
    #f.close()
