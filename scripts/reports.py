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

def _get_yield(demux_yields, sample):
    if demux_yields.shape[0] <= 0:
        return 0
    if demux_yields[demux_yields['Sample_ID'] == sample].shape[0] <= 0:
        return 0
    return demux_yields[demux_yields['Sample_ID'] == sample]['yield'].sum()

def get_status_data(samplesheets, config, prefix=None):
    if samplesheets.shape[0] <= 0:
        raise ValueError("No samples specified.")

    ACTION_to_PROGRAM = {
        'background': ['GATK', 'Platypus'],
        'trio': ['Varscan\ndenovo'],
        'tumornormal': ['Mutect', 'Varscan'],
        'demultiplex': ['bcl2fastq']}
    min_targets = 80

    if prefix is None:
        prefix = config['dirs']['prefix']

    # check if flowcell has been demuxed
    demux_yields = []
    check_flowcell_demuxed = dict()
    for run, g_run in samplesheets.groupby('run'):
        fp_yieldreport = join(prefix, config['dirs']['intermediate'], config['stepnames']['convert_illumina_report'], '%s.yield_report.pdf' % run)
        check_flowcell_demuxed[run] = exists(fp_yieldreport)
        if check_flowcell_demuxed[run]:
            # if run is already demultiplexed, obtain yield results
            demux_res = json.load(open('%s%s%s/%s/Stats/Stats.json' % (prefix, config['dirs']['intermediate'], config['stepnames']['demultiplex'], run), 'r'))
            for lane_res in demux_res['ConversionResults']:
                for sample_res in lane_res['DemuxResults']:
                    demux_yields.append({
                        'Lane': lane_res['LaneNumber'],
                        'run': demux_res['RunId'],
                        'Sample_ID': sample_res['SampleId'],
                        'yield': sample_res['Yield']})
    demux_yields = pd.DataFrame(demux_yields)

    results = []
    for project, g_project in samplesheets.groupby('Sample_Project'):
        for entity, g_entity in g_project.groupby('spike_entity_id'):
            if project in config['projects']:
                coverages = dict()
                for sample, g_sample in g_entity.groupby('Sample_ID'):
                    fp_coverage = join(prefix, config['dirs']['intermediate'], config['stepnames']['exome_coverage'], project, '%s.exome_coverage.csv' % sample)
                    if exists(fp_coverage):
                        coverage = pd.read_csv(fp_coverage, sep="\t")
                        coverages[sample] = coverage.loc[coverage['percent_cumulative'].apply(lambda x: abs(x-min_targets)).idxmin(), '#coverage']
                    else:
                        coverages[sample] = 0

                for action in config['projects'][project]['actions']:
                    fp_extracted = join(prefix, config['dirs']['intermediate'], config['stepnames']['snupy_extractsamples'], project, '%s.%s.extracted' % (entity, action))
                    status = exists(fp_extracted)
                    for program in ACTION_to_PROGRAM[action]:
                        for sample, g_sample in g_entity.groupby('Sample_ID'):
                            role = ""
                            if pd.notnull(g_sample['spike_entity_role'].unique()[0]):
                                role = g_sample['spike_entity_role'].unique()[0]
                            if role == 'sibling':
                                role = 'sibling: %s' % sample[len(entity):]
                            if ((action in ['trio']) and (role.split(':')[0] in ['patient', 'sibling'])) or\
                               ((action in ['tumornormal']) and (role.split(':')[0] in ['tumor'])) or\
                               ((action in ['background'])):
                                results.append({
                                    'project': project,
                                    'entity': entity,
                                    'sample': role,
                                    'action': action,
                                    'program': program,
                                    'status': int(status),
                                    'coverage': coverages[sample],
                                    'yield': _get_yield(demux_yields, sample)})
                    if action == 'trio':
                        expected_roles = set(['patient', 'father', 'mother'])
                    elif action == 'tumornormal':
                        expected_roles = set(['tumor', 'healthy'])
                    elif action == 'background':
                        expected_roles = set()
                    elif action == 'demultiplex':
                        expected_roles = set()
                    else:
                        raise ValueError("Unexpected action '%s'!" % action)
                    if len(expected_roles - set(g_entity['spike_entity_role'].unique())) > 0:
                        for role in (expected_roles - set(g_entity['spike_entity_role'].unique())):
                            results.append({
                                'project': project,
                                'entity': entity,
                                'sample': role,
                                'action': action,
                                'program': program,
                                'status': 0,
                                'coverage': -1,
                                'yield': _get_yield(demux_yields, role)})

                for sample, g_sample in g_entity.groupby('Sample_ID'):
                    role = ""
                    if pd.notnull(g_sample['spike_entity_role'].unique()[0]):
                        role = g_sample['spike_entity_role'].unique()[0]
                    if role == 'sibling':
                        role = 'sibling: %s' % sample[len(entity):]
                    results.append({
                        'project': project,
                        'entity': entity,
                        'sample': role,
                        'action': 'demultiplex',
                        'program': 'bcl2fastq',
                        'status': int(check_flowcell_demuxed[g_sample['run'].unique()[0]]),
                        'coverage': coverages[sample],
                        'yield': _get_yield(demux_yields, sample)})

        for sample, g_sample in g_project[pd.isnull(g_project['spike_entity_id'])].groupby('Sample_ID'):
            results.append({
                'project': project,
                'entity': "",
                'sample': sample,
                'action': 'demultiplex',
                'program': 'bcl2fastq',
                'status': int(check_flowcell_demuxed[g_sample['run'].unique()[0]]),
                'coverage': 0,
                'yield': _get_yield(demux_yields, sample)})
    #return results
    res = pd.pivot_table(
        data=pd.DataFrame(results),
        index=['project', 'entity', 'sample', 'yield', 'coverage'],
        columns=['action', 'program'],
        values='status',
        fill_value=-1,
    )

    # reorder actions
    res = res.loc[:, [
        ('demultiplex', 'bcl2fastq'),
        ('background', 'GATK'),
        ('background', 'Platypus'),
        ('tumornormal', 'Mutect'),
        ('tumornormal', 'Varscan'),
        ('trio', 'Varscan\ndenovo'),
    ]]

    return res


def write_status_update(filename, status_table, config, offset_rows=0, offset_cols=0, min_yield=5.0):
    """
    Parameters
    ----------
    filename : str
        Filepath to output Excel file.
    status_table : pd.DataFrame
        Result from function "get_status_data".
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
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet()

    # column header
    format_header = workbook.add_format({
        'rotation': 90,
        'bold': True,
        'valign': 'vcenter',
        'align': 'center'})
    for col, program in enumerate(['yield', 'coverage'] + list(status_table.columns.get_level_values(1))):
        worksheet.write(offset_rows, offset_cols+col+3, program.replace('yield', 'yield (GB)'), format_header)

    worksheet.set_column(offset_cols, offset_cols, 4)
    # set width of entity name
    worksheet.set_column(offset_cols+1, offset_cols+1, 12)
    worksheet.set_column(offset_cols+len(status_table.index.names)-2, offset_cols+len(status_table.index.names)+status_table.shape[1], 4)
    worksheet.set_row(offset_rows, 80)

    # row header: project
    format_project = workbook.add_format({
        'rotation': 90,
        'bold': True,
        'valign': 'vcenter',
        'align': 'center'})
    for pos, project in enumerate(status_table.index.levels[0]):
        start = list(status_table.index.labels[0]).index(pos)
        stop = len(status_table.index.labels[0]) - 1 - list(status_table.index.labels[0][::-1]).index(pos)
        col = chr(65+offset_cols)
        cellrange = '%s%i:%s%i' % (col, offset_rows+start+1+1, col, offset_rows+stop+1+1)
        worksheet.merge_range(cellrange, project.replace('_', '\n'), format_project)

    # row header: spike_entity_id
    format_spike_entity_id = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center'})
    for pos, spike_entity_id in enumerate(status_table.index.levels[1]):
        if spike_entity_id == "":
            continue
        start = list(status_table.index.labels[1]).index(pos)
        stop = len(status_table.index.labels[1]) - 1 - list(status_table.index.labels[1][::-1]).index(pos)
        col = chr(65+offset_cols+1)
        cellrange = '%s%i:%s%i' % (col, offset_rows+start+1+1, col, offset_rows+stop+1+1)
        worksheet.merge_range(cellrange, spike_entity_id, format_spike_entity_id)

    # row header: sample name / role
    format_spike_entity_role = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center'})
    format_spike_entity_role_missing = workbook.add_format({
        'valign': 'vcenter',
        'align': 'center',
        'font_color': '#ff0000'})
    for pos, role in enumerate(status_table.index.get_level_values(2)):
        if status_table.index.get_level_values(1)[pos] == "":
            cellrange = '%s%i:%s%i' % (chr(65+offset_cols+1), offset_rows+pos+1+1, chr(65+offset_cols+2), offset_rows+pos+1+1)
            worksheet.merge_range(cellrange, role, format_spike_entity_role)
        else:
            worksheet.write(offset_rows+1+pos, offset_cols+2, role, format_spike_entity_role if status_table.index[pos][-1] != -1 else format_spike_entity_role_missing)

    # yield
    format_yield = workbook.add_format({'align': 'center'})
    for row, yieldvalue in enumerate(status_table.index.get_level_values('yield')):
        worksheet.write(offset_rows+1+row, offset_cols+3, float('%.1f' % (yieldvalue / (1000**3))) if yieldvalue > 0 else '?', format_yield)
    format_yield_poor = workbook.add_format({'bg_color': '#ffcccc'})
    format_yield_good = workbook.add_format({'bg_color': '#ccffcc'})
    for pos, project in enumerate(status_table.index.levels[0]):
        start = list(status_table.index.labels[0]).index(pos)
        stop = len(status_table.index.labels[0]) - 1 - list(status_table.index.labels[0][::-1]).index(pos)
        #cellrange = '%s%i:%s%i' % (col, start+2, col, stop+2)
        cellrange = '%s%i:%s%i' % (chr(65+offset_cols+len(status_table.index.levels)-2), offset_rows+2+start,
                                   chr(65+offset_cols+len(status_table.index.levels)-2), offset_rows+2+stop)
        worksheet.conditional_format(cellrange, {'type': 'cell',
                                                 'criteria': '>=',
                                                 'value': 5.0,
                                                 'format': format_yield_good})
        worksheet.conditional_format(cellrange, {'type': 'cell',
                                                 'criteria': '<',
                                                 'value': 5.0,
                                                 'format': format_yield_poor})
    # coverage
    format_coverage = workbook.add_format({'align': 'center'})
    for row, coverage in enumerate(status_table.index.get_level_values('coverage')):
        worksheet.write(offset_rows+1+row, offset_cols+4, coverage if coverage > 0 else '?', format_coverage)
    format_coverage_poor = workbook.add_format({'bg_color': '#ffcccc'})
    format_coverage_good = workbook.add_format({'bg_color': '#ccffcc'})
    # conditional coloring based for coverage, based on project
    for pos, project in enumerate(status_table.index.levels[0]):
        start = list(status_table.index.labels[0]).index(pos)
        stop = len(status_table.index.labels[0]) - 1 - list(status_table.index.labels[0][::-1]).index(pos)
        #cellrange = '%s%i:%s%i' % (col, start+2, col, stop+2)
        cellrange = '%s%i:%s%i' % (chr(65+offset_cols+len(status_table.index.levels)-1), offset_rows+2+start,
                                   chr(65+offset_cols+len(status_table.index.levels)-1), offset_rows+2+stop)
        worksheet.conditional_format(cellrange, {'type': 'cell',
                                                 'criteria': '>=',
                                                 'value': get_min_coverage(project, config),
                                                 'format': format_coverage_good})
        worksheet.conditional_format(cellrange, {'type': 'cell',
                                                 'criteria': '<',
                                                 'value': get_min_coverage(project, config),
                                                 'format': format_coverage_poor})

    # status
    format_missing_sample = workbook.add_format({'align': 'center',
                                                 'font_color': '#ff0000'})
    for row, (idx_row, data) in enumerate(status_table.iterrows()):
        coverage = idx_row[-1]
        if coverage == -1:
            cellrange = '%s%i:%s%i' % (chr(65+offset_cols+4), row+offset_rows+2, chr(65+offset_cols+4+status_table.shape[1]), row+offset_rows+2)
            worksheet.merge_range(cellrange, 'missing sample', format_missing_sample)
        else:
            for col, (_, status) in enumerate(data.iteritems()):
                if pd.isnull(status) or (status == -1):
                    value = ""
                elif status == 1:
                    value = "done"
                elif status == 0:
                    value = "no"
                worksheet.write(offset_rows+1+row, offset_cols+5+col, value)
    format_processing_poor = workbook.add_format({'bg_color': '#ffcccc',
                                                  'font_color': '#ffcccc'})
    format_processing_good = workbook.add_format({'bg_color': '#ccffcc',
                                                  'font_color': '#ccffcc'})
    cellrange = '%s%i:%s%i' % (chr(65+len(idx_row)+offset_cols), offset_rows+2,
                               chr(65+len(idx_row)+status_table.shape[1]-1+offset_cols), status_table.shape[0]+offset_rows+2-1)
    worksheet.conditional_format(cellrange, {'type': 'cell',
                                             'criteria': '==',
                                             'value': '"done"',
                                             'format': format_processing_good})
    worksheet.conditional_format(cellrange, {'type': 'cell',
                                             'criteria': '==',
                                             'value': '"no"',
                                             'format': format_processing_poor})

    workbook.close()
