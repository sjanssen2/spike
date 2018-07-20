import socket
import os
from scripts.parse_samplesheet import get_fastq_filenames
from scripts.checks import check_illuminarun_complete

if socket.gethostname().startswith("hilbert") or socket.gethostname().startswith("murks"):
    # this loads a highly parallelized version of bcl2fastq compiled by HHU-HPC's guys
    shell.prefix("module load bcl2fastq;")

configfile: "config.yaml"
STEPNAME="05_demultiplex"


rule check_complete:
    input:
        fp_status=expand("{a}{b}/ImageAnalysis_Netcopy_complete.txt", a=config["dir_rawillumina"], b=config["run"])
    output:
        fp_check=expand("{dir}{stepname}-{run}-illuminarun_complete.txt", dir=config['dir_checks'], stepname=STEPNAME, run=config["run"])
    run:
        check_illuminarun_complete(input.fp_status[0], output.fp_check[0])


rule demultiplex:
    input:
        status=rules.check_complete.output
    output:
        # read expected fastq filenames from samplesheet and prefix with configured dir where to store them
        expand(map(lambda x: os.path.join(config["dir_intermediate"], STEPNAME, config["run"], x), get_fastq_filenames(os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])))),
    params:
        fp_samplesheet=os.path.join(config["dir_samplesheets"], "%s_ukd.csv" % config["run"])
    log:
        os.path.join(config["dir_logs"], config["run"], "%s.log" % STEPNAME)
    benchmark:
        os.path.join(config["dir_benchmarks"], config["run"], "%s.txt" % STEPNAME)
    threads: 100
    shell:
        "uname -a > {log} 2>&1; "
        "bcl2fastq"
        " --runfolder-dir {config[dir_rawillumina]}/{config[run]}/"
        " --output-dir {config[dir_intermediate]}/{STEPNAME}/{config[run]}/"
        " --ignore-missing-bcls"
        " --sample-sheet {params.fp_samplesheet}"
        " --loading-threads {threads}"
        " --processing-threads {threads}"
        " --writing-threads {threads}"
        " 2>> {log}"
