import os
from scripts.parse_samplesheet import get_lanes_for_sampleID
from glob import glob

configfile: "config.yaml"
STEPNAME="10_rejoin_samples"

rule rejoin_samples:
    input:
        "~/gpfs/Intermediate/10_rejoin_samples/180614_SN737_0438_BCC7MCACXX/AG_Remke/Chri_3/297_S1_R1.fastq.gz /gpfs/project/jansses/Intermediate/10_rejoin_samples/180614_SN737_0438_BCC7MCACXX/Alps/ALPS_66_a_S18_R1.fastq.gz"


rule test:
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

# snakemake -np ~/gpfs/Intermediate/10_rejoin_samples/180614_SN737_0438_BCC7MCACXX/AG_Remke/Chri_3/297_S1_R1.fastq.gz
# ~/gpfs/Intermediate/10_rejoin_samples/180614_SN737_0438_BCC7MCACXX/Alps/ALPS_66_a_S18_R{1,2}.fastq.gz

# rule rejoin_sample:
#     input:
#         glob("{dir}/Intermediate/05_demultiplex/{run}/{project}/{sample}_L*_R{direction}_001.fastq.gz")
#     output:
#         "{dir}/Intermediate/10_rejoin_samples/{run}/{project}/{sample}_R{direction,[1|2]}.fastq.gz"
#         #fastqs=expand(["%s_%s.fastq.gz" % (prefix, direction) for prefix in get_sample_fastqprefixes(os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])) for direction in ['R1', 'R2']])  # map(lambda x: os.path.join(config["dir_intermediate"], "05_demultiplex", config["run"], x), get_fastq_filenames(os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])))),
#         #fastqs=expand([os.path.join(config['dir_intermediate'], STEPNAME, config['run'], "%s_%s.fastq.gz" % (prefix, direction)) for prefix in get_sample_fastqprefixes(os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])) for direction in ['R1', 'R2']])  # map(lambda x: os.path.join(config["dir_intermediate"], "05_demultiplex", config["run"], x), get_fastq_filenames(os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])))),
#         #"kurt.txt"
#     shell:
#         "echo {input} > {output}"
#~/gpfs/Intermediate/10_rejoin_samples/180614_SN737_0438_BCC7MCACXX/AG_Remke/Chri_3/297_S1_R1.fastq.gz

# rule test:
#     input:
#         x=glob(expand('tes{word}_a*.txt', word=word))
#     output:
#         "tes{word}_b.txt"
#     shell:
#         "echo '{input.x}' > {output}"

# rule all:
#     input: dynamic("test_a{clusterid}.txt")
#
# rule plot:
#     input: "test_a{clusterid}.txt"
#     output: "test_b.txt"
#     shell:
#         "echo {input} > {output}"
