import os
from scripts.parse_samplesheet import get_sample_fastqprefixes
import glob
#from scripts.utils import load_modules
configfile: "config.yaml"

# DNA material of a sample might be split and loaded into several lanes to
# increase coverage. We here use information from SampleSheet and merge fastq.gz
# files from demultiplexing if necessary, otherwise we just use soft links

rule rejoin_samples:
    input:
        ['%s%s%s/%s/%s_%s.fastq.gz' % (config['dirs']['prefix'], config['dirs']['intermediate'], config['stepnames']['rejoin_samples'], config['run'], s, direction) for s in get_sample_fastqprefixes(os.path.join(
            config['dirs']['prefix'],
            config['dirs']['inputs'],
            config['dirs']['samplesheets'],
            "%s_ukd.csv" % config['run'])) for direction in config['directions']]


rule rejoin_sample:
    input:
        lambda wildcards: glob.glob('%s%s%s/%s/%s*%s_001.fastq.gz' % (wildcards.prefix, config["dirs"]["intermediate"], config["stepnames"]["demultiplex"], wildcards.run, wildcards.sample, wildcards.direction))
    benchmark:
        "{prefix}%s{run}/{sample}_R{direction}.benchmark" % config['dirs']['benchmarks']
    log:
        "{prefix}%s{run}/{sample}_R{direction}.log" % config['dirs']['logs']
    output:
        "{prefix}%s%s/{run,[^\/]+XX}/{sample, .*?}_{direction,R[1|2]}.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['rejoin_samples'])
    threads:
        1
    shell:
        'if [[ $(echo "{input}" | wc -w) -gt 1 ]]; then '
        # you can just concatenate multiple *.gz files into one, while
        # content when decompressed remains the same!
        'cat {input} > {output} 2> {log};'
        'else '
        'cp -l -v {input} {output} 2> {log}; '
        'chmod u+w {output} 2>> {log}; '
        'fi; '
