import socket
import os
from scripts.parse_samplesheet import get_fastq_filenames

if socket.gethostname().startswith("hpc-"):
    # this loads a highly parallelized version of bcl2fastq compiled by HHU-HPC's guys
    shell.prefix("module load bcl2fastq;")

configfile: "config.yaml"

rule demultiplex:
    output:
        # read expected fastq filenames from samplesheet and prefix with configured dir where to store them
        expand(map(lambda x: os.path.join(config["dir_demultiplexed"], x), get_fastq_filenames(os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"]))))
    params:
        fp_samplesheet=os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])
    shell:
        "bcl2fastq"
        " --runfolder-dir {config[dir_rawillumina]}/{config[run]}/"
        " --output-dir {config[dir_demultiplexed]}/{config[run]}/"
        " --ignore-missing-bcls"
        " --sample-sheet {params.fp_samplesheet}"
