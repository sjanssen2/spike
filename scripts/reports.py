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
import datetime
import getpass
import socket
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



ACTION_PROGRAMS = [
    {'action': 'background',
     'program': 'GATK',
     'fileending_snupy_extract': '.snp_indel.gatk',
     'fileending_spike_calls': '.gatk.snp_indel.vcf',
     'stepname_spike_calls': 'gatk_CombineVariants',
    },
    {'action': 'background',
     'program': 'Platypus',
     'fileending_snupy_extract': '.indel.ptp',
     'fileending_spike_calls': '.ptp.annotated.filtered.indels.vcf',
     'stepname_spike_calls': 'platypus_filtered',
    },
    {'action': 'tumornormal',
     'program': 'Varscan',
     'fileending_snupy_extract': '.somatic.varscan',
     'fileending_spike_calls': '.indel_snp.vcf',
     'stepname_spike_calls': 'merge_somatic',
    },
    {'action': 'tumornormal',
     'program': 'Mutect',
     'fileending_snupy_extract': '.somatic.mutect',
     'fileending_spike_calls': '.all_calls.vcf',
     'stepname_spike_calls': 'mutect',
    },
    {'action': 'trio',
     'program': 'Varscan\ndenovo',
     'fileending_snupy_extract': '.denovo.varscan',
     'fileending_spike_calls': '.var2denovo.vcf',
     'stepname_spike_calls': 'writing_headers',
    },
]


def _get_statusdata_demultiplex(samplesheets, prefix, config):
    demux_yields = []
    for flowcell in samplesheets['run'].unique():
        fp_yielddata = '%s%s%s/Data/%s.yield_data.csv' % (prefix, config['dirs']['intermediate'], config['stepnames']['yield_report'], flowcell)
        if exists(fp_yielddata):
            demux_yields.append(
                pd.read_csv(fp_yielddata, sep="\t").rename(columns={'Project': 'Sample_Project', 'Sample': 'Sample_ID', 'Yield': 'yield'})) #.set_index(['Project', 'Lane', 'Sample', 'Barcode sequence'])
    if len(demux_yields) <= 0:
        return pd.DataFrame()
    demux_yields = pd.concat(demux_yields, axis=0)
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
    if len(coverages) <= 0:
        return pd.DataFrame()
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
            for file_ending, action, program in [(ap['fileending_snupy_extract'], ap['action'], ap['program']) for ap in ACTION_PROGRAMS]:
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
    if len(results) <= 0:
        return pd.DataFrame()
    return pd.DataFrame(results).set_index(['Sample_Project', 'Sample_ID', 'action', 'program'])


def _get_statusdata_numberpassingcalls(samplesheets, prefix, config, RESULT_NOT_PRESENT):
    results = []
    for (sample_project, sample_id), meta in samplesheets.groupby(['Sample_Project', 'Sample_ID']):
        for file_ending, stepname, action, program in [(ap['fileending_spike_calls'], ap['stepname_spike_calls'], ap['action'], ap['program']) for ap in ACTION_PROGRAMS]:
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
    if len(results) <= 0:
        return pd.DataFrame()
    return pd.DataFrame(results).set_index(['Sample_Project', 'Sample_ID', 'action', 'program'])['number_calls']


def write_status_update(filename, samplesheets, config, prefix, offset_rows=0, offset_cols=0, min_yield=5.0, verbose=sys.stderr):
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

    print("Creating report", file=verbose)
    # obtain data
    print("1/5) gathering demuliplexing yields: ...", file=verbose, end="")
    data_yields = _get_statusdata_demultiplex(samplesheets, prefix, config)
    print(" done.\n2/5) gathering coverage: ...", file=verbose, end="")
    data_coverage = _get_statusdata_coverage(samplesheets, prefix, config)
    print(" done.\n3/5) gathering snupy extraction status: ...", file=verbose, end="")
    data_snupy = _get_statusdata_snupyextracted(samplesheets, prefix, config)
    print(" done.\n4/5) gathering number of PASSing calls: ...", file=verbose, end="")
    data_calls = _get_statusdata_numberpassingcalls(samplesheets, prefix, config, RESULT_NOT_PRESENT)
    print("done.\n5/5) generating Excel output: ...", file=verbose, end="")

    # start creating the Excel result
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet()

    format_good = workbook.add_format({'bg_color': '#ccffcc'})
    format_bad = workbook.add_format({'bg_color': '#ffcccc'})

    # date information
    format_info = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center',
        'font_size': 9})
    cellrange = '%s%i:%s%i' % (chr(65+offset_cols), offset_rows+1, chr(65+offset_cols+3), offset_rows+1)
    info_username = getpass.getuser()
    info_now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    info_machine = socket.gethostname()
    worksheet.merge_range(cellrange, ('status report created\nat %s\nby %s\non %s' % (info_now, info_username, info_machine)),format_info)

    # header
    format_header = workbook.add_format({
        'rotation': 90,
        'bold': True,
        'valign': 'vcenter',
        'align': 'center'})
    worksheet.set_row(offset_rows, 80)
    for i, caption in enumerate(['yield (MB)', 'coverage'] + [ap['program'] for ap in ACTION_PROGRAMS]):
        worksheet.write(offset_rows, offset_cols+4+i, caption, format_header)
    format_spike_seqdate = workbook.add_format({
        'align': 'center',
        'valign': 'vcenter',
        'font_size': 8})
    worksheet.write(offset_rows, offset_cols+6+len(ACTION_PROGRAMS), 'sequenced at', format_spike_seqdate)

    worksheet.freeze_panes(offset_rows+1, offset_cols+4)

    # body
    format_project = workbook.add_format({
        'rotation': 90,
        'bold': True,
        'valign': 'vcenter',
        'align': 'center'})
    format_spike_entity_id = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center'})
    format_spike_sampleID = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center'})
    format_spike_entity_role_missing = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center',
        'font_color': '#ff0000'})
    row = offset_rows+1
    for sample_project, grp_project in samplesheets.groupby('Sample_Project'):
        cellrange = '%s%i:%s%i' % (chr(65+offset_cols), row+1, chr(65+offset_cols), row+len(grp_project.groupby(['spike_entity_id', 'Sample_ID'])))
        worksheet.merge_range(cellrange, sample_project.replace('_', '\n'), format_project)
        worksheet.set_column(offset_cols, offset_cols, 4)

        # add in lines to indicate missing samples, e.g. for trios that are incomplete
        missing_samples = []
        if (sample_project in config['projects']) and ('actions' in config['projects'][sample_project]):
            if 'trio' in config['projects'][sample_project]['actions']:
                for spike_entity_id, grp_spike_entity_group in grp_project.groupby('spike_entity_id'):
                    for role in ['patient', 'mother', 'father']:
                        if grp_spike_entity_group[grp_spike_entity_group['spike_entity_role'] == role].shape[0] <= 0:
                            missing_samples.append({
                                'spike_entity_id': spike_entity_id,
                                'Sample_ID': role,
                                'spike_entity_role': role,
                                'missing': True,
                            })

        # groupby excludes NaNs, thus I have to hack: replace NaN by "" here and
        # reset to np.nan within the loop
        for spike_entity_group, grp_spike_entity_group in pd.concat([grp_project, pd.DataFrame(missing_samples)], sort=False).fillna(value={'spike_entity_id': ''}).groupby('spike_entity_id'):
            cellrange = '%s%i:%s%i' % (chr(65+offset_cols+1), row+1, chr(65+offset_cols+1), row+len(grp_spike_entity_group.groupby('Sample_ID')))
            if spike_entity_group != "":
                if len(set(cellrange.split(':'))) > 1:
                    worksheet.merge_range(cellrange, spike_entity_group, format_spike_entity_id)
                else:
                    worksheet.write(cellrange.split(':')[0], spike_entity_group, format_spike_entity_id)
            else:
                spike_entity_group = np.nan
            worksheet.set_column(offset_cols+1, offset_cols+1, 10)

            for nr_sample_id, (sample_id, grp_sample_id) in enumerate(grp_spike_entity_group.sort_values(by='spike_entity_role').groupby('Sample_ID')):
                worksheet.set_column(offset_cols+2, offset_cols+2, 4)
                role = grp_sample_id['spike_entity_role'].iloc[0]
                is_missing = ('missing' in grp_sample_id.columns) and (grp_sample_id[grp_sample_id['missing'] == np.True_].shape[0] > 0)

                # sample_ID, extend field if no spike_entity_group or spike_entity_role is given
                col_start = offset_cols+2
                col_end = offset_cols+2
                if pd.isnull(spike_entity_group):
                    col_start -= 1
                if pd.isnull(role):
                    col_end += 1
                cellrange = '%s%i:%s%i' % (chr(65+col_start), row+1, chr(65+col_end), row+1)
                sample_id_value = sample_id
                # if sample_id starts with name of the entity group, we are using "..." to make it visually more pleasing
                if pd.notnull(spike_entity_group) and sample_id_value.startswith(spike_entity_group):
                    sample_id_value = '%s' % sample_id[len(spike_entity_group):]
                frmt = format_spike_sampleID
                if is_missing:
                    sample_id_value = '?'
                    frmt = format_spike_entity_role_missing
                if len(set(cellrange.split(':'))) > 1:
                    worksheet.merge_range(cellrange, sample_id_value, frmt)
                else:
                    worksheet.write(cellrange.split(':')[0], sample_id_value, frmt)

                # spike_entity_role
                if pd.notnull(role):
                    worksheet.write(row, offset_cols+3, str(role), frmt)

                if is_missing:
                    cellrange = '%s%i:%s%i' % (chr(65+offset_cols+4), row+1, chr(65+offset_cols+3+2+len(ACTION_PROGRAMS)+1), row+1)
                    worksheet.merge_range(cellrange, "missing sample", format_spike_entity_role_missing)
                else:
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
                    if pd.notnull(role):
                        frmt = format_bad
                        value_coverage = "missing"
                        if (sample_project, sample_id) in data_coverage.index:
                            value_coverage = data_coverage.loc[sample_project, sample_id]
                            if value_coverage >= get_min_coverage(sample_project, config):
                                frmt = format_good
                        worksheet.write(row, offset_cols+5, value_coverage, frmt)
                        worksheet.set_column(offset_cols+5, offset_cols+5, 4)

                    for i, (action, program) in enumerate([(ap['action'], ap['program']) for ap in ACTION_PROGRAMS]):
                        value_numcalls = ""
                        frmt = None
                        if (sample_project, sample_id, action, program) in data_calls:
                            value_numcalls = data_calls.loc[sample_project, sample_id, action, program]
                        if (sample_project, sample_id, action, program) in data_snupy.index:
                            if data_snupy.loc[sample_project, sample_id, action, program]['status']:
                                frmt = format_good
                            else:
                                frmt = format_bad
                        if value_numcalls == RESULT_NOT_PRESENT:
                            value_numcalls = 'vcf missing'
                            frmt = format_bad
                        if frmt is not None:
                            worksheet.write(row, offset_cols+6+i, value_numcalls, frmt)

                    # sequencing date
                    worksheet.write(row, offset_cols+6+len(ACTION_PROGRAMS), ' / '.join(sorted(map(
                        lambda x: datetime.datetime.strptime('20%s' % x.split('_')[0], '%Y%m%d').strftime("%Y-%m-%d"), grp_sample_id['run'].unique()))), format_spike_seqdate)
                    worksheet.set_column(offset_cols+6+len(ACTION_PROGRAMS), offset_cols+6+len(ACTION_PROGRAMS), 16)

                row += 1
    print("done.\n", file=verbose, end="")
    workbook.close()


def collect_yield_data(dir_flowcell):
    """Composes yield report for (potentially) split demulitplexing.

    Notes
    -----
    People used different length barcodes in the same lane / flowcell.
    This should in general be avoided, due to potential clashes of barcodes.
    Therefore Illumina's bcl2fastq fails to process sample sheets formatted
    like this. I overcome this issue by splitting the samplesheet and running
    bcl2fast multiple times independently. However, this requires complicated
    logic especially for stats for the undetermined reads.

    Parameters
    ----------
    dir_flowcell : str
        Filepath to demultiplexing results in split fashion.

    Returns
    -------
    3-tuple of pd.DataFrame : (Flowcell Summary, Lane Summary, Top Unknown Barcodes)
    Formatting is required before data will resemble Illumina's original yield report.
    """
    lane_summary = []
    lane_meta = []
    cluster_meta = []
    unknown_barcodes = []
    for dir_part in glob(join(dir_flowcell, 'part_*')):
        clusters = []
        for fp_fastqsummary in glob(join(dir_part, 'Stats/FastqSummaryF*L*.txt')):
            cluster = pd.read_csv(fp_fastqsummary, sep="\t")
            perc_pf_clusters = pd.concat([cluster.groupby(['SampleNumber'])['NumberOfReadsPF'].sum(),
                                          cluster.groupby(['SampleNumber'])['NumberOfReadsRaw'].sum()], axis=1)
            perc_pf_clusters['Lane'] = fp_fastqsummary.split('/')[-1].split('.')[0].split('L')[-1]
            clusters.append(perc_pf_clusters)
        clusters = pd.concat(clusters).reset_index().set_index(['Lane', 'SampleNumber'])
        cluster_meta.append(clusters)

        meta_samples = pd.read_csv(join(dir_part, 'Stats/AdapterTrimming.txt'), sep="\t", usecols=[0,2,3,4])
        meta_samples = meta_samples.iloc[:meta_samples[meta_samples['Lane'].apply(lambda x: x.startswith('Lane:'))].index.max()-2,:].drop_duplicates()
        meta_samples.set_index(['Lane', 'Sample Id'], inplace=True)

        fp_json = join(dir_part, 'Stats/Stats.json')
        part_stats = json.load(open(fp_json, 'r'))
        # sample numbers in FastqSummary files match S-idx numbers, which are only increased if sample is not seen before, independent on lane
        sample_numbers = dict()
        for res_conv in part_stats['ConversionResults']:
            numq30bases = 0
            sumQualityScore = 0
            for res_demux in res_conv['DemuxResults']:
                q30bases = sum([res_metrics['YieldQ30'] for res_metrics in res_demux['ReadMetrics']])
                qualityScore = sum([res_metrics['QualityScoreSum'] for res_metrics in res_demux['ReadMetrics']])
                if res_demux['SampleId'] not in sample_numbers:
                    sample_numbers[res_demux['SampleId']] = len(sample_numbers)+1
                sample_number = sample_numbers[res_demux['SampleId']]
                lane_summary.append({
                    'Lane': res_conv['LaneNumber'],
                    'Project': meta_samples.loc[str(res_conv['LaneNumber']), res_demux['SampleId']]['Project'],
                    'Sample': res_demux['SampleId'],
                    'Barcode sequence': res_demux['IndexMetrics'][0]['IndexSequence'],
                    'Barcode length': len(res_demux['IndexMetrics'][0]['IndexSequence']),
                    'PF Clusters': res_demux['NumberReads'],
                    '% of the lane': res_demux['NumberReads'] / res_conv['TotalClustersPF'],
                    '% Perfect barcode': sum([res_idx['MismatchCounts']['0'] for res_idx in res_demux['IndexMetrics']]) / res_demux['NumberReads'] if res_demux['NumberReads'] > 0 else 0,
                    '% One mismatch barcode': sum([res_idx['MismatchCounts']['1'] for res_idx in res_demux['IndexMetrics']]) / res_demux['NumberReads'] if res_demux['NumberReads'] > 0 else 0,
                    'Yield': res_demux['Yield'],
                    '% PF Clusters': clusters.loc[str(res_conv['LaneNumber']), sample_number]['NumberOfReadsPF'] / clusters.loc[str(res_conv['LaneNumber']), sample_number]['NumberOfReadsRaw'] if clusters.loc[str(res_conv['LaneNumber']), sample_number]['NumberOfReadsRaw'] > 0 else 0,
                    '% >= Q30 bases': q30bases / res_demux['Yield'] if res_demux['Yield'] > 0 else 0,
                    'Q30 bases': q30bases,
                    'QualityScoreSum': qualityScore,
                    'Mean Quality Score': sum([res_metrics['QualityScoreSum'] for res_metrics in res_demux['ReadMetrics']]) / res_demux['Yield'] if res_demux['Yield'] > 0 else 0,
                    'Sample_Number': sample_number,
                })
                numq30bases += q30bases
                sumQualityScore += qualityScore
            numq30bases += sum([res_metrics['YieldQ30'] for res_metrics in res_conv['Undetermined']['ReadMetrics']])
            sumQualityScore += sum([res_metrics['QualityScoreSum'] for res_metrics in res_conv['Undetermined']['ReadMetrics']])

            lane_meta.append({
                'Lane': res_conv['LaneNumber'],
                "TotalClustersRaw" : res_conv['TotalClustersRaw'],
                "TotalClustersPF" : res_conv['TotalClustersPF'],
                "Yield": res_conv['Yield'],
                "YieldQ30": numq30bases,
                "QualityScoreSum": sumQualityScore,
                "Flowcell": part_stats['Flowcell'],
            })
        for res_unknown in part_stats['UnknownBarcodes']:
            for barcode in res_unknown['Barcodes'].keys():
                unknown_barcodes.append({
                    'Lane': res_unknown['Lane'],
                    'Barcode sequence': barcode,
                    'Count': res_unknown['Barcodes'][barcode],
                    'Barcode length': len(barcode)})

    lane_meta = pd.DataFrame(lane_meta).drop_duplicates().set_index('Lane')
    lane_meta['run'] = dir_flowcell.split('/')[-1]
    lane_summary = pd.DataFrame(lane_summary)
    cluster_meta = pd.concat(cluster_meta)

    undetermined = []
    for lane, clst_lane in lane_summary.groupby('Lane'):
        undetermined.append({
            'Lane': lane,
            'Project': 'default',
            'Sample': 'Undetermined',
            'Barcode sequence': 'unknown',
            'PF Clusters': lane_meta.loc[lane, 'TotalClustersPF'] - clst_lane['PF Clusters'].sum(),
            '% of the lane': (lane_meta.loc[lane, 'TotalClustersPF'] - clst_lane['PF Clusters'].sum()) / lane_meta.loc[lane, 'TotalClustersPF'],
            '% Perfect barcode': 1,
            '% One mismatch barcode': np.nan,
            'Yield': lane_meta.loc[lane, 'Yield'] - clst_lane['Yield'].sum(),
            '% PF Clusters': (lane_meta.loc[lane, 'TotalClustersPF'] - cluster_meta.loc[str(lane), range(1, cluster_meta.index.get_level_values('SampleNumber').max()), :]['NumberOfReadsPF'].sum()) / (lane_meta.loc[lane, 'TotalClustersRaw'] - cluster_meta.loc[str(lane), range(1, cluster_meta.index.get_level_values('SampleNumber').max()), :]['NumberOfReadsRaw'].sum()),
            '% >= Q30 bases': (lane_meta.loc[lane, 'YieldQ30'] - lane_summary[lane_summary['Lane'] == lane]['Q30 bases'].sum()) / (lane_meta.loc[lane, 'Yield'] - lane_summary[lane_summary['Lane'] == lane]['Yield'].sum()),
            'Mean Quality Score': (lane_meta.loc[lane, 'QualityScoreSum'] - lane_summary[lane_summary['Lane'] == lane]['QualityScoreSum'].sum()) / (lane_meta.loc[lane, 'Yield'] - lane_summary[lane_summary['Lane'] == lane]['Yield'].sum()),
        })

    lane_summary = pd.concat([lane_summary, pd.DataFrame(undetermined)], sort=False)

    # traverse unknown barcodes and filter those that are actually used by samples
    # this is non-trivial, because due to splitting, used barcodes might have different sizes!
    unknown_barcodes = pd.DataFrame(unknown_barcodes)
    idx_remove = []
    for lane in unknown_barcodes['Lane'].astype(int).unique():
        lane_known_bcs = lane_summary[(lane_summary['Lane'] == lane) & (lane_summary['Barcode sequence'] != 'unknown')][['Barcode sequence', 'Barcode length']]
        lane_unknown_bcs = unknown_barcodes[unknown_barcodes['Lane'] == lane]
        remove = set()
        for len_known, known in lane_known_bcs.groupby('Barcode length'):
            for len_unknown, unknown in lane_unknown_bcs.groupby('Barcode length'):
                if len_known == len_unknown:
                    remove |= set(unknown['Barcode sequence']) & set(known['Barcode sequence'])
                elif len_known < len_unknown:
                    for _, bc_unknown in unknown['Barcode sequence'].iteritems():
                        if bc_unknown[:int(len_known)] in set(known['Barcode sequence']):
                            remove |= set([bc_unknown])
                elif len_known > len_unknown:
                    remove |= set(unknown['Barcode sequence']) & set(known['Barcode sequence'].apply(lambda x: x[:int(len_unknown)]))
        idx_remove.extend(lane_unknown_bcs[lane_unknown_bcs['Barcode sequence'].isin(remove)].index)
    top_unknown_barcodes = unknown_barcodes.loc[set(unknown_barcodes.index) - set(idx_remove),:]

    return lane_meta, lane_summary, top_unknown_barcodes


def create_html_yield_report(fp_yield_report, lane_meta, lane_summary, top_unknown_barcodes, config):
    """Creates a HTML yield report.

    Parameters
    ----------
    fp_yield_report : str
        Filepath to resulting HTML file.
    lane_meta : pd.DataFrame
        Result of collect_yield_data, first component:
        Flowcell summary statistics, basically # clusters.
    lane_summary : pd.DataFrame
        Result of collect_yield_data, second component:
        Sample demultiplexing information.
    top_unknown_barcodes : pd.DataFrame
        Result of collect_yield_data, third component:
        Infos about high abundant unused barcodes.
    config : dict
        Snakemakes configuration dictionary.
    """
    #dir_flowcell : str
    #    Filepath to demultiplexing results in split fashion.

    out = '<html>\n<head>\n'
    #out += '<link rel="stylesheet" href="Report.css" type="text/css">\n'
    out += '<style>\n'
    out += 'body {font-size: 100%; font-family:monospace;}\n'
    out += 'table#ReportTable {border-width: 1px 1px 1px 1px; border-collapse: collapse;}\n'
    out += 'table#ReportTable td {text-align:right; padding: 0.3em;}\n'
    out += 'thead { display: table-header-group }\n'
    out += 'tfoot { display: table-row-group }\n'
    out += 'tr { page-break-inside: avoid }\n'
    out += '</style>\n'
    out += '<title>%s</title>\n' % lane_meta['run'].unique()[0]
    out += '</head>\n<body>'
    out += "%s / [all projects] / [all samples] / [all barcodes]" % lane_meta['Flowcell'].unique()[0]


    out += "<h2>Flowcell Summary</h2>"
    fc_summary = lane_meta.sum()[['TotalClustersRaw', 'TotalClustersPF', 'Yield']].to_frame().T
    fc_summary.rename(columns={'TotalClustersRaw': 'Clusters (Raw)',
                               'TotalClustersPF': 'Clusters(PF)',
                               'Yield': 'Yield (MBases)'}, inplace=True)
    for col in ['Clusters (Raw)', 'Clusters(PF)']:
        fc_summary[col] = fc_summary[col].apply(lambda x: '{:,}'.format(x))
    fc_summary['Yield (MBases)'] = fc_summary['Yield (MBases)'].apply(lambda x: '{:,}'.format(int(x/1000000)))
    out += fc_summary.to_html(index=False, table_id='ReportTable', justify='center')


    out += "<h2>Lane Summary</h2>"
    x = lane_summary.sort_values(['Lane', 'Sample_Number'])[['Lane', 'Project', 'Sample', 'Barcode sequence', 'PF Clusters', '% of the lane', '% Perfect barcode', '% One mismatch barcode', 'Yield', '% PF Clusters', '% >= Q30 bases', 'Mean Quality Score']].rename(columns={'Yield': 'Yield (Mbases)'})
    x['PF Clusters'] = x['PF Clusters'].apply(lambda x: '{:,}'.format(x))
    for col in ['% of the lane', '% Perfect barcode', '% PF Clusters', '% >= Q30 bases', '% One mismatch barcode']:
        x[col] = x[col].apply(lambda x: '%.2f' % (x*100))
    x['Yield (Mbases)'] = x['Yield (Mbases)'].apply(lambda x: '{:,}'.format(int(x/1000000)))
    x['Mean Quality Score'] = x['Mean Quality Score'].apply(lambda x: '%.2f' % x)
    x['% One mismatch barcode'] = x['% One mismatch barcode'].replace('nan', 'NaN')
    out += x.to_html(index=False, table_id='ReportTable', justify='center')


    out += "<h2>Top Unknown Barcodes</h2>"
    topX = 10
    lanes = []
    for lane, lane_barcodes in top_unknown_barcodes.groupby('Lane'):
        x = lane_barcodes.sort_values('Count', ascending=False).iloc[:topX][['Lane', 'Count', 'Barcode sequence']].rename(columns={'Barcode sequence': 'Sequence'})#.reset_index().set_index(['Lane', 'Count'])
        x.index = range(1,topX+1)[:x.shape[0]]
        x['Count'] = x['Count'].apply(lambda x: '{:,}'.format(x))
        lanes.append(x)
    topunknown = pd.concat(lanes, axis=1)
    out += '<table border="1" class="dataframe" id="ReportTable">\n<thead>\n<tr style="text-align: center;">\n'
    for col in topunknown.columns:
        out += '<th>%s</th>\n' % col
    out += '</tr>\n</thead>\n<tbody>\n'
    for i, row in topunknown.iterrows():
        out += '<tr>\n'
        for col, value in row.iteritems():
            if col == 'Lane':
                if i == 1:
                    out += '<th rowspan=%s>%s</th>\n' % (topX, value)
            else:
                out += '<td>%s</td>\n' % value
        out += '</tr>\n'
    out += '</tbody>\n</table>\n'


    out += "Report generated by %s.</body>\n</html>" % config['name_program']

    with open(fp_yield_report, 'w') as f:
        f.write(out)
