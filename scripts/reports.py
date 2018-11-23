from os.path import join, exists, dirname, basename
from os import makedirs
import sys
import pandas as pd
from glob import glob
import seaborn as sns
import numpy as np
from scipy import stats
import xlsxwriter
import matplotlib.pyplot as plt
from scripts.parse_samplesheet import get_min_coverage
import json
import requests
from requests.auth import HTTPBasicAuth
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
plt.switch_backend('Agg')


def report_undertermined_filesizes(fp_filesizes, fp_output, fp_error,
                                   zscorethreshold=1):
    # read all data
    fps_sizes = glob(join(dirname(fp_filesizes), '*.txt'))
    pds_sizes = []
    for fp_size in fps_sizes:
        data = pd.read_csv(
            fp_size, sep="\t", names=["filesize", "filename", "status"],
            index_col=1)
        # mark given read as isme=True while all other data in the dir
        # are isme=False
        data['isme'] = fp_filesizes in fp_size
        data['filesize'] /= 1024**3
        pds_sizes.append(data)
    pd_sizes = pd.concat(pds_sizes)

    # compute z-score against non-bad known runs
    pd_sizes['z-score'] = np.nan
    idx_nonbad = pd_sizes[pd_sizes['status'] != 'bad'].index
    pd_sizes.loc[idx_nonbad, 'z-score'] = stats.zscore(
        pd_sizes.loc[idx_nonbad, 'filesize'])

    # plot figure
    fig = plt.figure()
    ax = sns.distplot(
        pd_sizes[(pd_sizes['isme'] == np.False_) &
                 (pd_sizes['status'] != 'bad')]['filesize'],
        kde=False, rug=False, color="black", label='known runs')
    ax = sns.distplot(
        pd_sizes[(pd_sizes['isme'] == np.False_) &
                 (pd_sizes['status'] == 'bad')]['filesize'],
        kde=False, rug=False, color="red", label='bad runs')
    ax = sns.distplot(
        pd_sizes[pd_sizes['isme'] == np.True_]['filesize'],
        kde=False, rug=True, color="green", label='this run')
    _ = ax.set_ylabel('number of files')
    _ = ax.set_xlabel('file-size in GB')
    ax.set_title('run %s' % basename(fp_filesizes)[:-4])
    ax.legend()

    # raise error if current run contains surprisingly large undetermined
    # filesize
    if pd_sizes[(pd_sizes['isme'] == np.True_) &
                (pd_sizes['status'] == 'unknown')]['z-score'].max() > zscorethreshold:
        ax.set_title('ERROR: %s' % ax.get_title())
        fig.savefig(fp_error, bbox_inches='tight')
        raise ValueError(
            ("Compared to known historic runs, your run contains surprisingly "
             "(z-score > %f) large file(s) of undetermined reads. You will find"
             " an supporting image at '%s'. Please do the following things:\n"
             "1. discuss with lab personal about the quality of the run.\n"
             "2. should you decide to keep going with this run, mark file "
             "status (3rd column) in file '%s' as 'good'.\n"
             "3. for future automatic considerations, mark file status (3rd "
             "column) as 'bad' if you have decided to abort processing due to"
             " too low quality (z-score kind of averages about known values)."
             ) % (zscorethreshold, fp_error, fp_filesizes))
    else:
        fig.savefig(fp_output, bbox_inches='tight')


def report_exome_coverage(
    fps_sample, fp_plot,
    min_coverage=30, min_targets=80, coverage_cutoff=200):
    """Creates an exome coverage plot for multiple samples.

    Parameters
    ----------
    fps_sample : [str]
        A list of file-paths with coverage data in csv format.
    fp_plot : str
        Filepath of output graph.
    min_coverage : int
        Default: 30.
        An arbitraty threshold of minimal coverage that we expect.
        A vertical dashed line is drawn at this value.
    min_targets : float
        Default: 80.
        An arbitraty threshold of minimal targets that we expect to be covered.
        A horizontal dashed line is drawn at this value.
    coverage_cutoff : float
        Default: 200.
        Rightmost coverage cut-off value where X-axis is limited.

    Raises
    ------
    ValueError : If one of the sample's coverage falls below expected
    thresholds.
    """
    # Usually we aim for a 30X coverage on 80% of the sites.

    fig, ax = plt.subplots()
    ax.axhline(y=min_targets, xmin=0, xmax=coverage_cutoff, color='gray',
               linestyle='--')
    ax.axvline(x=min_coverage, ymin=0, ymax=100, color='gray', linestyle='--')

    samples_below_coverage_threshold = []
    for fp_sample in fps_sample:
        coverage = pd.read_csv(fp_sample, sep="\t")
        samplename = fp_sample.split('/')[-1].split('.')[0]
        linewidth = 1
        if coverage[coverage['#coverage'] == min_coverage]['percent_cumulative'].min() < min_targets:
            linewidth = 4
            samples_below_coverage_threshold.append(samplename)
        ax.plot(coverage['#coverage'],
                coverage['percent_cumulative'],
                label=samplename,
                linewidth=linewidth)

    ax.set_xlim((0, coverage_cutoff))
    ax.set_xlabel('Read Coverage')
    ax.set_ylabel('Targeted Exome Bases')
    ax.legend()

    if len(samples_below_coverage_threshold) > 0:
        fp_plot = fp_plot.replace('.pdf', '.error.pdf')

    fig.savefig(fp_plot, bbox_inches='tight')

    if len(samples_below_coverage_threshold) > 0:
        raise ValueError(
            "The following %i sample(s) have coverage below expected "
            "thresholds. Please discuss with project PIs on how to proceed. "
            "Maybe, samples need to be re-sequenced.\n\t%s\nYou will find more"
            " information in the generated coverage plot '%s'." % (
                len(samples_below_coverage_threshold),
                '\n\t'.join(samples_below_coverage_threshold),
                fp_plot))


def _get_statusdata_demultiplex(samplesheets, prefix, config):
    demux_yields = []
    for flowcell in samplesheets['run'].unique():
        fp_demuxstats = '%s%s%s/%s/Stats/Stats.json' % (prefix, config['dirs']['intermediate'], config['stepnames']['demultiplex'], flowcell)
        if exists(fp_demuxstats):
            demux_res = json.load(open(fp_demuxstats, 'r'))
            for lane_res in demux_res['ConversionResults']:
                for sample_res in lane_res['DemuxResults']:
                    inferred_sample_project = samplesheets[
                        (samplesheets['Sample_ID'] == sample_res['SampleId']) &
                        (samplesheets['Lane'] == lane_res['LaneNumber']) &
                        (samplesheets['run'] == demux_res['RunId'])]['Sample_Project'].unique()
                    if len(inferred_sample_project) > 1:
                        raise ValueError('Conflicting Sample_Project description for sample %s in run %s.' % (sample_res['SampleId'], demux_res['RunId']))
                    if len(inferred_sample_project) < 1:
                        # this means, that samples in the demux stats but not in the provided samplesheets are
                        # NOT part of the result
                        continue
                    demux_yields.append({
                        'Lane': lane_res['LaneNumber'],
                        'run': demux_res['RunId'],
                        'Sample_Project': inferred_sample_project[0],
                        'Sample_ID': sample_res['SampleId'],
                        'yield': sample_res['Yield']})
    return pd.DataFrame(demux_yields).groupby(['Sample_Project', 'Sample_ID'])['yield'].sum()


def _get_statusdata_coverage(samplesheets, prefix, config, min_targets=80):
    coverages = []
    for (sample_project, sample_id), meta in samplesheets.groupby(['Sample_Project', 'Sample_ID']):
        fp_coverage = join(prefix, config['dirs']['intermediate'], config['stepnames']['exome_coverage'], sample_project, '%s.exome_coverage.csv' % sample_id)
        if exists(fp_coverage):
            coverage = pd.read_csv(fp_coverage, sep="\t")
            coverages.append({
                'Sample_Project': sample_project,
                'Sample_ID': sample_id,
                'coverage': coverage.loc[coverage['percent_cumulative'].apply(lambda x: abs(x-min_targets)).idxmin(), '#coverage']})
    return pd.DataFrame(coverages).set_index(['Sample_Project', 'Sample_ID'])['coverage']


def _get_statusdata_snupyextracted(samplesheets, prefix, config):
    results = []
    for sample_project, meta in samplesheets.groupby('Sample_Project'):
        if (sample_project not in config['projects']) or ('snupy' not in config['projects'][sample_project]) or ('project_id' not in config['projects'][sample_project]['snupy']):
            # project in config file is not properly configure for snupy!
            continue
        r = requests.get('%s/experiments/%s.json' % (config['credentials']['snupy']['host'], config['projects'][sample_project]['snupy']['project_id']),
            auth=HTTPBasicAuth(config['credentials']['snupy']['username'], config['credentials']['snupy']['password']),
            verify=False)
        assert(r.headers.get('status') == '200 OK')
        samples = [sample['name'] for sample in r.json()['samples']]

        for sample_id, meta_sample in meta.groupby('Sample_ID'):
            for file_ending, action, program in [('.snp_indel.gatk', 'background', 'GATK'),
                                                 ('.indel.ptp', 'background', 'Platypus'),
                                                 ('.denovo.varscan', 'trio', 'Varscan\ndenovo'),
                                                 ('.somatic.varscan', 'tumornormal', 'Varscan'),
                                                 ('.somatic.mutect', 'tumornormal', 'Mutect')]:
                # in some cases "sample name" hold spike_entity_id, in others Sample_ID
                entity = sample_id
                runs = '+'.join(sorted(meta_sample['run'].unique()))
                if (action == 'trio'):
                    if meta_sample['spike_entity_role'].unique()[0] == 'patient':
                        entity = meta_sample['spike_entity_id'].iloc[0]
                        runs = '+'.join(sorted(samplesheets[samplesheets['spike_entity_id'] == meta_sample['spike_entity_id'].iloc[0]]['run'].unique()))
                if (action == 'tumornormal'):
                    if meta_sample['spike_entity_role'].unique()[0] == 'tumor':
                        entity = meta_sample['spike_entity_id'].iloc[0]
                        runs = '+'.join(sorted(samplesheets[samplesheets['spike_entity_id'] == meta_sample['spike_entity_id'].iloc[0]]['run'].unique()))
                name = '%s_%s/%s%s' % (runs, sample_project, entity, file_ending)

                if action in config['projects'][sample_project]['actions']:
                    if ((action == 'trio') and (meta_sample['spike_entity_role'].iloc[0] in ['patient', 'sibling'])) or\
                       ((action == 'background')) or\
                       ((action == 'tumornormal') and (meta_sample['spike_entity_role'].iloc[0].startswith('tumor'))):
                        results.append({
                            'Sample_Project': sample_project,
                            'Sample_ID': sample_id,
                            'action': action,
                            'program': program,
                            'status': name in samples,
                            'snupy_sample_name': name
                        })

    return pd.DataFrame(results).set_index(['Sample_Project', 'Sample_ID', 'action', 'program'])


def _get_statusdata_numberpassingcalls(samplesheets, prefix, config, RESULT_NOT_PRESENT):
    results = []
    for (sample_project, sample_id), meta in samplesheets.groupby(['Sample_Project', 'Sample_ID']):
        for file_ending, stepname, action, program in [
            ('.gatk.snp_indel.vcf', 'gatk_CombineVariants', 'background', 'GATK'),
            ('.ptp.annotated.filtered.indels.vcf', 'platypus_filtered', 'background', 'Platypus'),
            ('.var2denovo.vcf', 'writing_headers', 'trio', 'Varscan\ndenovo'),
            ('.indel_snp.vcf', 'merge_somatic', 'tumornormal', 'Varscan'),
            ('.all_calls.vcf', 'mutect', 'tumornormal', 'Mutect')
            ]:
            name = sample_id
            if (action == 'trio'):
                if meta['spike_entity_role'].unique()[0] == 'patient':
                    name = meta['spike_entity_id'].iloc[0]
            if (action == 'tumornormal'):
                if meta['spike_entity_role'].unique()[0] == 'tumor':
                    name = meta['spike_entity_id'].iloc[0]
            fp_vcf = '%s%s%s/%s/%s%s' % (prefix, config['dirs']['intermediate'], config['stepnames'][stepname], sample_project, name, file_ending)

            nr_calls = RESULT_NOT_PRESENT
            if exists(fp_vcf):
                nr_calls = pd.read_csv(fp_vcf, comment='#', sep="\t", dtype=str, header=None, usecols=[6], squeeze=True).value_counts()['PASS']

            if (sample_project in config['projects']) and (action in config['projects'][sample_project]['actions']):
                if ((action == 'trio') and (meta['spike_entity_role'].iloc[0] in ['patient', 'sibling'])) or\
                   ((action == 'background')) or\
                   ((action == 'tumornormal') and (meta['spike_entity_role'].iloc[0].startswith('tumor'))):
                    results.append({
                        'Sample_Project': sample_project,
                        'Sample_ID': sample_id,
                        'action': action,
                        'program': program,
                        'number_calls': nr_calls,
                    })
    return pd.DataFrame(results).set_index(['Sample_Project', 'Sample_ID', 'action', 'program'])['number_calls']


def write_status_update_v2(filename, samplesheets, config, prefix, offset_rows=0, offset_cols=0, min_yield=5.0):
    """
    Parameters
    ----------
    filename : str
        Filepath to output Excel file.
    samplesheets : pd.DataFrame
        The global samplesheets.
    config : dict()
        Snakemake configuration object.
    offset_rows : int
        Default: 0
        Number if rows to leave blank on top.
    offset_cols : int
        Default: 0
        Number if columns to leave blank on the left.
    min_yield : float
        Default: 5.0
        Threshold when to color yield falling below this value in red.
        Note: I don't know what a good default looks like :-/
    """
    RESULT_NOT_PRESENT = -5

    # obtain data
    print("1/5) gathering demuliplexing yields: ...", file=sys.stderr, end="")
    data_yields = _get_statusdata_demultiplex(samplesheets, prefix, config)
    print(" done.\n2/5) gathering coverage: ...", file=sys.stderr, end="")
    data_coverage = _get_statusdata_coverage(samplesheets, prefix, config)
    print(" done.\n3/5) gathering snupy extraction status: ...", file=sys.stderr, end="")
    data_snupy = _get_statusdata_snupyextracted(samplesheets, prefix, config)
    print(" done.\n4/5) gathering number of PASSing calls: ...", file=sys.stderr, end="")
    data_calls = _get_statusdata_numberpassingcalls(samplesheets, prefix, config, RESULT_NOT_PRESENT)
    print("done.\n5/5) generating Excel output: ...", file=sys.stderr, end="")

    # start creating the Excel result
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet()

    format_good = workbook.add_format({'bg_color': '#ccffcc'})
    format_bad = workbook.add_format({'bg_color': '#ffcccc'})
    actionsprograms = [
        ('background', 'GATK'),
        ('background', 'Platypus'),
        ('tumornormal', 'Varscan'),
        ('tumornormal', 'Mutect'),
        ('trio', 'Varscan\ndenovo')]

    # header
    format_header = workbook.add_format({
        'rotation': 90,
        'bold': True,
        'valign': 'vcenter',
        'align': 'center'})
    worksheet.set_row(offset_rows, 80)
    for i, caption in enumerate(['yield (MB)', 'coverage'] + list(map(lambda x: x[-1], actionsprograms))):
        worksheet.write(offset_rows, offset_cols+4+i, caption, format_header)

    # body
    format_project = workbook.add_format({
        'rotation': 90,
        'bold': True,
        'valign': 'vcenter',
        'align': 'center'})
    format_spike_entity_id = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center'})
    row = offset_rows+1
    for sample_project, grp_project in samplesheets.groupby('Sample_Project'):
        cellrange = '%s%i:%s%i' % (chr(65+offset_cols), row+1, chr(65+offset_cols), row+len(grp_project.groupby(['spike_entity_id', 'Sample_ID'])))
        worksheet.merge_range(cellrange, sample_project.replace('_', '\n'), format_project)
        worksheet.set_column(offset_cols, offset_cols, 4)

        for spike_entity_group, grp_spike_entity_group in grp_project.fillna(value={'spike_entity_id': 'NONE'}).groupby('spike_entity_id'):
            #worksheet.write(row, offset_cols+1, spike_entity_group)
            cellrange = '%s%i:%s%i' % (chr(65+offset_cols+1), row+1, chr(65+offset_cols+1), row+len(grp_spike_entity_group.groupby('Sample_ID')))
            if (spike_entity_group != "NONE") and (len(set(cellrange.split(':'))) > 1):
                worksheet.merge_range(cellrange, spike_entity_group, format_spike_entity_id)
            else:
                worksheet.set_column(offset_cols+1, offset_cols+1, 12)

            for sample_id, grp_sample_id in grp_spike_entity_group.sort_values(by='spike_entity_role').groupby('Sample_ID'):
                if spike_entity_group != "NONE":
                    worksheet.write(row, offset_cols+2, sample_id, format_spike_entity_id)
                else:
                    cellrange = '%s%i:%s%i' % (chr(65+offset_cols+1), row+1, chr(65+offset_cols+2), row+1)
                    worksheet.merge_range(cellrange, sample_id, format_spike_entity_id)
                worksheet.set_column(offset_cols+2, offset_cols+2, 12)

                worksheet.write(row, offset_cols+3, str(grp_sample_id['spike_entity_role'].fillna('').iloc[0]), format_spike_entity_id)

                # demultiplexing yield
                frmt = format_bad
                value_yield = "missing"
                if (sample_project, sample_id) in data_yields.index:
                    value_yield = float('%.1f' % (int(data_yields.loc[sample_project, sample_id]) / (1000**3)))
                    if value_yield >= 5.0:
                        frmt = format_good
                worksheet.write(row, offset_cols+4, value_yield, frmt)
                worksheet.set_column(offset_cols+4, offset_cols+4, 4)

                # coverage
                frmt = format_bad
                value_coverage = "missing"
                if (sample_project, sample_id) in data_coverage.index:
                    value_coverage = data_coverage.loc[sample_project, sample_id]
                    if value_coverage >= get_min_coverage(sample_project, config):
                        frmt = format_good
                worksheet.write(row, offset_cols+5, value_coverage, frmt)
                worksheet.set_column(offset_cols+5, offset_cols+5, 4)

                for i, (action, program) in enumerate(actionsprograms):
                    value_numcalls = ""
                    frmt = None
                    if (sample_project, sample_id, action, program) in data_calls:
                        value_numcalls = data_calls.loc[sample_project, sample_id, action, program]
                    if (sample_project, sample_id, action, program) in data_snupy.index:
                        if data_snupy.loc[sample_project, sample_id, action, program]['status']:
                            frmt = format_good
                        else:
                            frmt = format_bad
                    if frmt is not None:
                        worksheet.write(row, offset_cols+6+i, value_numcalls, frmt)

                row += 1
    print("done.\n", file=sys.stderr, end="")
    workbook.close()
