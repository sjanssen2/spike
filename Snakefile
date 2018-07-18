import socket
import os
from scripts.parse_samplesheet import get_fastq_filenames

if socket.gethostname().startswith("hilbert") or socket.gethostname().startswith("hpc-"):
    # this loads a highly parallelized version of bcl2fastq compiled by HHU-HPC's guys
    shell.prefix("module load bcl2fastq;")

configfile: "config.yaml"

rule demultiplex:
    output:
        # read expected fastq filenames from samplesheet and prefix with configured dir where to store them
        expand(map(lambda x: os.path.join(config["dir_demultiplexed"], config["run"], x), get_fastq_filenames(os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"]))))
    params:
        fp_samplesheet=os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])
    log:
        os.path.join(config["dir_logs"], config["run"], "05_demultiplexed.log")
    benchmark:
        os.path.join(config["dir_benchmarks"], config["run"], "05_demultiplexed.txt")
    shell:
        "uname -a > {log} 2>&1; "
        "bcl2fastq"
        " --runfolder-dir {config[dir_rawillumina]}/{config[run]}/"
        " --output-dir {config[dir_demultiplexed]}/{config[run]}/"
        " --ignore-missing-bcls"
        " --sample-sheet {params.fp_samplesheet}"
        " 2>> {log}"
