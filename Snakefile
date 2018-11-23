import os
import socket
import glob
import pandas as pd

from scripts.parse_samplesheet import parse_samplesheet, get_role, get_reference_genome, get_reference_knowns, get_reference_exometrack, get_species, get_reference_varscan_somatic, get_global_samplesheets
from scripts.parse_samplesheet import get_trios, get_tumorNormalPairs, get_samples, get_bwa_mem_header, get_demux_samples, get_projects_with_exomecoverage, get_rejoin_input, get_xenograft_hybridreference, get_xenograft_stepname, get_persamplefastq_samples, get_min_coverage
from scripts.utils import exclude_sample
from scripts.reports import report_undertermined_filesizes, report_exome_coverage, write_status_update
from scripts.convert_platypus import annotate
from scripts.snupy import upload_to_snupy, extractsamples


if socket.gethostname().startswith("hilbert") or socket.gethostname().startswith("murks"):
    # I now prefer to use the bcl2fastq version shipped via conda!
    # this loads a highly parallelized version of bcl2fastq compiled by HHU-HPC's guys
    # shell.prefix("module load bcl2fastq;")
    shell.prefix("module load R;")

configfile: "config.yaml"
SAMPLESHEETS = get_global_samplesheets(os.path.join(config['dirs']['prefix'], config['dirs']['inputs'], config['dirs']['samplesheets']))

include: "rules/demultiplex/Snakefile"
include: "rules/backup/Snakefile"
include: "rules/rejoin_samples/Snakefile"
include: "rules/trim/Snakefile"
include: "rules/xenograft/Snakefile"
include: "rules/map/Snakefile"
include: "rules/gatk/Snakefile"
include: "rules/platypus/Snakefile"
include: "rules/varscan/Snakefile"
include: "rules/freec/Snakefile"
include: "rules/mutect/Snakefile"
include: "rules/snupy/Snakefile"
include: "rules/excavator/Snakefile"


localrules: check_complete, aggregate_undetermined_filesizes, check_undetermined_filesizes, convert_illumina_report, check_coverage, xenograft_check, correct_genotypes_somatic, varscan_fpfilter_somatic, somatic_FPfilter, vcf_annotate, merge_somatic_mus_musculus, merge_somatic_homo_sapiens, writing_headers, merge_vcfs, varscan_filter_INDEL, varscan_processSomatic

rule all:
    input:
        # create a yield report per run as one of the checkpoints for the wetlab crew
        # don't create yield report for per sample fastq projects!
        yield_report=["%s%s%s/%s.yield_report.pdf" % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['convert_illumina_report'], run)
                      for run in set(SAMPLESHEETS[pd.notnull(SAMPLESHEETS['Lane'])]['run'].unique())],

        # create backup for each run
        # don't backup per sample fastq runs
        # backup=["%s%s%s.%s.done" % (config['dirs']['prefix'], config['dirs']['checks'], run, config['stepnames']['backup_validate'])
        #         for run in set(SAMPLESHEETS['run'].unique()) - set(get_persamplefastq_samples(SAMPLESHEETS, config))],

        # demultiplex all samples for projects that ONLY need to demultiplex, e.g. AG_Remke
        demux=['%s%s%s/%s' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['demultiplex'], run)
               for run in get_demux_samples(SAMPLESHEETS, config)],

        # compute exome coverage for exome samples
        exome_coverage=["%s%s%s/%s/%s.exome_coverage.csv" % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['exome_coverage'], project, sample)
                        for _, (project, sample) in SAMPLESHEETS[SAMPLESHEETS['Sample_Project'].isin(get_projects_with_exomecoverage(config))][['Sample_Project', 'Sample_ID']].iterrows()],

        # upload and extract called variants to Snupy
        background=['%s%s%s/%s.background.extracted' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['snupy_extractsamples'], entity)
                    for entity in set(map(lambda x: '%s/%s' % (x['Sample_Project'], x['spike_entity_id']), get_samples(SAMPLESHEETS, config)))],
        tumornormal=['%s%s%s/%s/%s.tumornormal.extracted' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['snupy_extractsamples'], pair['Sample_Project'], pair['spike_entity_id']) for pair in get_tumorNormalPairs(SAMPLESHEETS, config)],
        trio_calling=['%s%s%s/%s/%s.trio.extracted' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['snupy_extractsamples'], trio['Sample_Project'], trio['spike_entity_id']) for trio in get_trios(SAMPLESHEETS, config)],

        # tumornormal calling for complete tumor/normal pairs of all runs
        # tumornormal_freec=['%s%s%s/%s/%s/tumor.pileup.freec_BAF.txt'  % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['freec'],        pair['Sample_Project'], pair['spike_entity_id']) for pair in get_tumorNormalPairs(SAMPLESHEETS, config)],

rule status:
    output:
        "status_update.xlsx"
    run:
        write_status_update(output[0], SAMPLESHEETS, config, prefix=config['dirs']['prefix'])

# onerror:
#     print("An error occurred")
#     shell("mail -s 'an error occurred' stefan.m.janssen@gmail.com < {log}")
