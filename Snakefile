import os
import socket
import glob
import pandas as pd

from scripts.parse_samplesheet import parse_samplesheet, get_role, get_reference_genome, get_reference_knowns, get_reference_exometrack, get_species, get_reference_varscan_somatic, get_global_samplesheets, split_samplesheets
from scripts.parse_samplesheet import get_trios, get_tumorNormalPairs, get_samples, get_bwa_mem_header, get_demux_samples, get_rejoin_input, get_xenograft_hybridreference, get_xenograft_stepname, get_min_coverage, get_genepanels, get_reverse_file, get_kind_of_run
from scripts.utils import exclude_sample, sample_to_biom, merge_samples
from scripts.reports import report_undertermined_filesizes, report_exome_coverage, get_status_data, write_status_update, create_html_yield_report, collect_yield_data, get_gene_panel_coverage
from scripts.convert_platypus import annotate
from scripts.snupy import upload_to_snupy, extractsamples


if socket.gethostname().startswith("hilbert") or socket.gethostname().startswith("murks"):
    # I now prefer to use the bcl2fastq version shipped via conda!
    # this loads a highly parallelized version of bcl2fastq compiled by HHU-HPC's guys
    # shell.prefix("module load bcl2fastq;")
    shell.prefix("module load R;")
    shell.prefix("module load kraken/2.0.7;")

SNUPY_INSTANCE = 'hhu'
configfile: "config.yaml"
SAMPLESHEETS = get_global_samplesheets(os.path.join(config['dirs']['prefix'], config['dirs']['inputs'], config['dirs']['samplesheets']), config)

if False:
    # subset for Keimbahn project:
    # last two conditions are to include samples of other projects, which also are used for specific roles in the Keimbahn project, see config.yaml for details.
    SAMPLESHEETS = SAMPLESHEETS[(SAMPLESHEETS['Sample_Project'] == 'Keimbahn') |
                                (SAMPLESHEETS['Sample_ID'].isin(['HL_ini', 'HL_rem']) & (SAMPLESHEETS['Sample_Project'] == 'ALL_Study1_Hauer')) |
                                (SAMPLESHEETS['Sample_ID'].isin(['ALPS_60', 'ALPS_60a', 'ALPS_60b']) & (SAMPLESHEETS['Sample_Project'] == 'Alps'))]

print("%i samples in %i projects." % (SAMPLESHEETS['Sample_ID'].unique().shape[0], SAMPLESHEETS['Sample_Project'].unique().shape[0]), file=sys.stderr)

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
include: "rules/analyses/Snakefile"


localrules: check_complete, aggregate_undetermined_filesizes, check_undetermined_filesizes, convert_illumina_report, check_coverage, xenograft_check, correct_genotypes_somatic, varscan_fpfilter_somatic, somatic_FPfilter, vcf_annotate, merge_somatic_mus_musculus, merge_somatic_homo_sapiens, writing_headers, merge_vcfs, varscan_filter_INDEL, varscan_processSomatic, split_demultiplex, yield_report

rule all:
    input:
        # demultiplex all samples for projects that ONLY need to demultiplex, e.g. AG_Remke and create yield reports
        demux=[res for run in get_demux_samples(SAMPLESHEETS, config) for res in [
               '%s%s%s/%s' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['join_demultiplex'], run),
               '%s%s%s/%s.yield_report.pdf' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['convert_illumina_report'], run)]],

        # STATISTICS ON BAM FILES of "PrintReads"
        # compute exome coverage for exome samples
        exome_coverage=["%s%s%s/%s.exome_coverage.csv" % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['exome_coverage'], sample)
                        for sample in set(map(lambda x: x['sample'], get_samples(SAMPLESHEETS, config)))],
        # compute bamstat exome samples
        bamstats=["%s%s%s/%s.bamstat.tsv" % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['bamstat'], sample)
                  for sample in set(map(lambda x: x['sample'], get_samples(SAMPLESHEETS, config)))],
        # compute gene panel coverage
        genepanels=get_genepanels(SAMPLESHEETS, config, config['dirs']['prefix']),

        # upload and extract called variants to Snupy
        background=['%s%s%s/%s/%s.background.extracted' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['snupy_extractsamples'], SNUPY_INSTANCE, entity)
                    for entity in set(map(lambda x: '%s/%s' % (x['Sample_Project'], x['spike_entity_id']), get_samples(SAMPLESHEETS, config)))],
        tumornormal=['%s%s%s/%s/%s/%s.tumornormal.extracted' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['snupy_extractsamples'], SNUPY_INSTANCE, pair['Sample_Project'], pair['spike_entity_id'])
                     for pair in get_tumorNormalPairs(SAMPLESHEETS, config)],
        trio_calling=['%s%s%s/%s/%s/%s.trio.extracted' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['snupy_extractsamples'], SNUPY_INSTANCE, trio['Sample_Project'], trio['spike_entity_id'])
                      for trio in get_trios(SAMPLESHEETS, config)],

        # STATISTICS ON BACKGROUND GATK
        biom_background=['%s%s%s/%s.gatk.snp_indel.biom' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['biom_gatkbackground'], sample)
                         for sample in set(map(lambda x: x['sample'], get_samples(SAMPLESHEETS, config)))],

rule backup:
    input:
        # create backup for each run
        # don't backup per sample fastq runs
        backup=["%s%s%s/%s.done" % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['backup_validate'], run)
                for run in set(SAMPLESHEETS['run'].unique()) - set(SAMPLESHEETS[pd.isnull(SAMPLESHEETS['Lane'])]['run'].unique())],

rule trim_16s:
    input:
        ['%s%s%s/%s_R%i.fastq.gz' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['remove_16s_primer'], row['fastq-prefix'], direction)
         for idx, row in SAMPLESHEETS.iterrows()
         for direction in [1, 2]]

rule status:
    output:
        "status_update.xlsx"
    run:
        status_data = get_status_data(SAMPLESHEETS, config, snupy_instance=SNUPY_INSTANCE)
        write_status_update(status_data, output[0], SAMPLESHEETS, config)

# onerror:
#     print("An error occurred")
#     shell("mail -s 'an error occurred' stefan.m.janssen@gmail.com < {log}")
