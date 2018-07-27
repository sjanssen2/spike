import os
import socket
import glob

from scripts.parse_samplesheet import get_sample_fastqprefixes
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
#
# rule all_rejoin_samples:
#     input:
#         ['%s%s%s/%s/%s_%s.fastq.gz' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['rejoin_samples'], config['run'], s, direction) for s in get_sample_fastqprefixes(os.path.join(
#             config['dirs']['prefix'],
#             config['dirs']['inputs'],
#             config['dirs']['samplesheets'],
#             "%s_ukd.csv" % config['run'])) for direction in config['directions']]
