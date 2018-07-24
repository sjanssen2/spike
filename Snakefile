import os
from scripts.parse_samplesheet import get_sample_fastqprefixes
import glob
#from scripts.utils import load_modules
configfile: "config.yaml"

# DNA material of a sample might be split and loaded into several lanes to
# increase coverage. We here use information from SampleSheet and merge fastq.gz
# files from demultiplexing if necessary, otherwise we just use soft links

rule all_trim:
    input:
        [ '%s%s%s/%s/%s/%s_%s.fastq.gz' % (
            config['dirs']['prefix'],
            config['dirs']['intermediate'],
            config['stepnames']['trim'],
            config['run'],
            pair,
            sample,
            direction) for sample in get_sample_fastqprefixes('~/gpfs/Inputs/SampleSheets/180614_SN737_0438_BCC7MCACXX_ukd.csv') for direction in ["R1", "R2"] for pair in ['Paired', 'Unpaired']]

rule all_rejoin_samples:
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

rule trim:
    input:
        forward="{prefix}%s%s/{run,[^\/]+XX}/{sample, .*?}_R1.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['rejoin_samples']),
        reverse="{prefix}%s%s/{run,[^\/]+XX}/{sample, .*?}_R2.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['rejoin_samples'])
    output:
        pairedforward="{prefix}%s%s/{run,[^\/]+XX}/Paired/{sample, .*?}_R1.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['trim']),
        unpairedforward="{prefix}%s%s/{run,[^\/]+XX}/Unpaired/{sample, .*?}_R1.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['trim']),
        pairedreverse="{prefix}%s%s/{run,[^\/]+XX}/Paired/{sample, .*?}_R2.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['trim']),
        unpairedreverse="{prefix}%s%s/{run,[^\/]+XX}/Unpaired/{sample, .*?}_R2.fastq.gz" % (config['dirs']['intermediate'], config['stepnames']['trim'])
    log:
        "{prefix}%s{run}/{sample}_trimmomatic.log" % config['dirs']['logs']
    benchmark:
        "{prefix}%s{run}/{sample}_trimmomatic.benchmark" % config['dirs']['logs']
    conda:
      "envs/spike_trim.yaml"
    threads:
        16
    shell:
        "java"
        " -Xmx4g"
        " -XX:ParallelGCThreads={threads}"
        " -jar ${{CONDA_PREFIX}}/share/trimmomatic-0.33-1/trimmomatic.jar"
        " PE -threads {threads} -phred33"
        " -trimlog {log}.trimlog"
        " {input.forward}"
        " {input.reverse}"
        " {output.pairedforward}"
        " {output.unpairedforward}"
        " {output.pairedreverse}"
        " {output.unpairedreverse}"
        " ILLUMINACLIP:${{CONDA_PREFIX}}/share/trimmomatic-0.33-1/adapters/TruSeq3-PE.fa:2:30:10 CROP:99 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36"
        " > {log} 2>&1"
