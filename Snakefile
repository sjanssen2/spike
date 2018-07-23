import os
from scripts.parse_samplesheet import get_lanes_for_sampleID, get_sample_fastqprefixes
from glob import glob

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
        lambda wildcards: ["%s%s%s/%s/%s/%s%s_S%s_L%03i_R%s_001.fastq.gz" % (
            wildcards.prefix,
            config['dirs']['intermediate'],
            config['stepnames']['demultiplex'],
            wildcards.run,
            wildcards.project,
            wildcards.samplegrp,
            wildcards.sample,
            wildcards.sidx,
            lane,
            wildcards.direction) for lane in get_lanes_for_sampleID(
                os.path.join(
                    wildcards.prefix,
                    config["dirs"]["inputs"],
                    config["dirs"]["samplesheets"],
                    "%s_ukd.csv" % config["run"]),
                wildcards.samplegrp, wildcards.sample, wildcards.sidx)]
    benchmark:
        "{prefix}%s{run}/{project}_`echo '{samplegrp}' | tr '/' '_'`{sample}_S{sidx}_R{direction}.benchmark" % config['dirs']['benchmarks']
    log:
        "{prefix}%s{run}/{project}_`echo '{samplegrp}' | tr '/' '_'`{sample}_S{sidx}_R{direction}.log" % config['dirs']['logs']
    output:
        "{prefix}%s%s/{run,[^/]*}/{project,[^/]*}/{samplegrp,.*?/{0,1}}{sample,[^/]+}_S{sidx,\d*}_R{direction,[1|2]}.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['rejoin_samples'])
    shell:
        "if [[ `echo '{input}' | wc -w` -gt 1 ]]; then "
        "zcat {input} | gzip > {output} 2> {log};"
        "else "
        "ln -v -s {input} {output} > {log} 2>&1; "
        "fi; "
