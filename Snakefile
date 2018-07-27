import os
import socket
import glob

from scripts.parse_samplesheet import get_sample_fastqprefixes, parse_samplesheet
#from scripts.utils import load_modules

if socket.gethostname().startswith("hilbert") or socket.gethostname().startswith("murks"):
    # this loads a highly parallelized version of bcl2fastq compiled by HHU-HPC's guys
    shell.prefix("module load bcl2fastq;")

configfile: "config.yaml"
include: "rules/demultiplex/Snakefile"
include: "rules/rejoin_samples/Snakefile"
include: "rules/trim/Snakefile"
include: "rules/map/Snakefile"
include: "rules/gatk/Snakefile"
include: "rules/platypus/Snakefile"

def _get_samples(fp_samplesheet):
    ss = parse_samplesheet(fp_samplesheet)
    samples = {'%s/%s_S%s' % (row['Sample_Project'], row['Sample_ID'], row['s-idx']) for idx, row in ss.iterrows()}
    return [sample for sample in samples if sample.startswith('Alps/ALPS_66')]

rule all:
    input:
        # the yield report as one of the checkpoints for the wetlab crew
        yield_report="%s%s%s/%s.yield_report.pdf" % (config['dirs']['prefix'], config['dirs']['reports'], config['run'], config['run']),

        # snv calling against background, this is the "main pipeline"
        snvs_background_ptp=['%s%s%s/%s/%s.ptp.annotated.filtered.indels.vcf' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['platypus_filtered'], config['run'], sample)
                             for sample in _get_samples("%s%s%s%s_ukd.csv" % (config['dirs']['prefix'], config['dirs']['inputs'], config['dirs']['samplesheets'], config['run']))],
        snvs_background_gatk=['%s%s%s/%s/%s.gatk.%ssnp_indel.vcf' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['gatk_CombineVariants'], config['run'], sample, isrelax)
                              for sample in _get_samples("%s%s%s%s_ukd.csv" % (config['dirs']['prefix'], config['dirs']['inputs'], config['dirs']['samplesheets'], config['run']))
                              for isrelax in ['', 'relax.']]

        

rule all_trim:
    input:
        [ '%s%s%s/%s/%s/%s_%s.fastq.gz' % (
            config['dirs']['prefix'],
            config['dirs']['intermediate'],
            config['stepnames']['trim'],
            config['run'],
            pair,
            sample,
            direction) for sample in get_sample_fastqprefixes('%s%s%s%s_ukd.csv' % (config['dirs']['prefix'], config['dirs']['inputs'], config['dirs']['samplesheets'], config['run'])) for direction in ["R1", "R2"] for pair in ['Paired', 'Unpaired']]

# rule all_rejoin_samples:
#     input:
#         ['%s%s%s/%s/%s_%s.fastq.gz' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['rejoin_samples'], config['run'], s, direction) for s in get_sample_fastqprefixes(os.path.join(
#             config['dirs']['prefix'],
#             config['dirs']['inputs'],
#             config['dirs']['samplesheets'],
#             "%s_ukd.csv" % config['run'])) for direction in config['directions']]
