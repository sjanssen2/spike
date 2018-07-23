import os
from scripts.parse_samplesheet import get_sample_fastqprefixes, get_laneSplitInputs
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
        lambda wildcards: get_laneSplitInputs(
            wildcards,
            os.path.join(config["dirs"]["inputs"], config["dirs"]["samplesheets"]),
            os.path.join(config["dirs"]["intermediate"], config["stepnames"]["demultiplex"]))
    benchmark:
        "{prefix}%s{run}/{sample}_R{direction}.benchmark" % config['dirs']['benchmarks']
    # log:
    #     "{prefix}%s{run}/{sample}_R{direction}.log" % config['dirs']['logs']
    output:
        "{prefix}%s%s/{run,[^\/]+}/{sample}_{direction,R[1|2]}.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['rejoin_samples'])
    threads:
        1
    shell:
        'if [[ $(echo "{input}" | wc -w) -gt 1 ]]; then '
        # you can just concatenate multiple *.gz files into one, while
        # content when decompressed remains the same!
        'cat {input} > {output};'
        'else '
        'cp -l -v {input} {output}; '
        'chmod u+w {output}; '
        'fi; '
