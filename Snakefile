import socket
import os
from scripts.parse_samplesheet import get_fastq_filenames

if socket.gethostname().startswith("hpc-"):
    # this loads a highly parallelized version of bcl2fastq compiled by HHU-HPC's guys
    shell.prefix("module load bcl2fastq;")

configfile: "config.yaml"

rule demultiplex:
    input:
        s=os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])
    output:
        "test.txt"
    script:
        print(input().__dict__)
    #shell:
    #    "echo 'Hallo Welt' > test.txt; bcl2fastq --version >> test.txt 2>&1"
