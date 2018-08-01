import os
import socket
import glob

from scripts.parse_samplesheet import get_sample_fastqprefixes, parse_samplesheet, get_sample_names, get_xenograft_host
from scripts.utils import exclude_sample


if socket.gethostname().startswith("hilbert") or socket.gethostname().startswith("murks"):
    # this loads a highly parallelized version of bcl2fastq compiled by HHU-HPC's guys
    shell.prefix("module load bcl2fastq;")

configfile: "config.yaml"
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

EXCLUDE_SAMPLES = ['Maus_Hauer', 'Fischer']

rule all:
    input:
        # the yield report as one of the checkpoints for the wetlab crew
        yield_report="%s%s%s/%s.yield_report.pdf" % (config['dirs']['prefix'], config['dirs']['reports'], config['run'], config['run']),

        # snv calling against background, this is the "main pipeline"
        snvs_background_ptp=['%s%s%s/%s/%s.ptp.annotated.filtered.indels.vcf' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['platypus_filtered'], config['run'], sample)
                             for sample in get_sample_names("%s%s%s%s_ukd.csv" % (config['dirs']['prefix'], config['dirs']['inputs'], config['dirs']['samplesheets'], config['run']))
                             if not exclude_sample(sample, EXCLUDE_SAMPLES)],
        snvs_background_gatk=['%s%s%s/%s/%s.gatk.%ssnp_indel.vcf' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['gatk_CombineVariants'], config['run'], sample, isrelax)
                              for sample in get_sample_names("%s%s%s%s_ukd.csv" % (config['dirs']['prefix'], config['dirs']['inputs'], config['dirs']['samplesheets'], config['run']))
                              for isrelax in ['', 'relax.']
                              if not exclude_sample(sample, EXCLUDE_SAMPLES)],

        coverage_report='%s%s%s/%s.exome_coverage.pdf' % (config['dirs']['prefix'], config['dirs']['reports'], config['run'], config['run']),

        trio=['%s%s%s/%s/%s.var2denovo.vcf' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['writing_headers'], 'Alps', 'ALPS_66')],

        somatic_freec=['%s%s%s/%s/%s/tumor.pileup.gz_BAF.txt' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['freec'], 'Fischer_Geron', 'hum_leuk_unknown')],
        somatic_mutect=['%s%s%s/%s/%s.all_calls.csv' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['mutect'], 'Fischer_Geron', 'hum_leuk_unknown')],

        backup="%s%s%s.%s.done" % (config['dirs']['prefix'], config['dirs']['checks'], config['run'], config['stepnames']['backup_validate']),

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
