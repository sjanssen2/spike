## Requirements
You might have to install (mini|ana)conda https://docs.conda.io/en/latest/miniconda.html, create an environment https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-environments and install snakemake https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda first.

# Read this first
`Snakemake` greatly supports execution of a pipeline on a grid compute system with very little overhead. However, it will help reading Snakemake's own documentation to get familiar with with basic concepts and commands: https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution

We are here listing our best practice experience running `spike` on the `High Performance Computing` (HPC) system at University of Duesseldorf: https://www.zim.hhu.de/high-performance-computing.html Great people support users on the HPC via https://rocketchat.hhu.de/ or by phone or email and it is definitively worth reading their Wiki https://wiki.hhu.de/display/HPC/Wissenschaftliches+Hochleistungs-Rechnen+am+ZIM **before** submitting thousands of jobs to their grid.

# Compose the snakemake master command
Here, I am giving some background information about the multiple flags you should set when executing a `spike` pipeline via snakemake.

## `--cluster-config cluster.json`
Resources on shared computer clusters are not infinite and the job of the admins in to ensure that multiple users get their fair share. This is only possible, if each user behaves nicely, i.e. spends some efforts to find small upper boundaries for their programs needs in terms of memory, CPUs/cores and execution time. Naturally, you are greedy, but admins counteract by putting greedy jobs late in the queue, i.e. you have to wait longer until your greedy job is executed compared to a small job. More details are in Snakemakes documentation. https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html?highlight=cluster-status#cluster-configuration

`spike` has ~70 different rules, each rule running a different command with different resource demands, depending on typical input sizes. Thus, finding working upper bounds required some experience, which is now distilled in the file https://github.com/sjanssen2/spike/blob/master/cluster.json First entry is the default, which will be applied to every rule if no further matching definition can be found. 
 - In this default rule, we define the `account` to which the used service units get "billed" to. I once had to contact HPC admins and request to create this new account for spike.
 - The HPC has several queues into which your jobs can be submitted. For `spike`, we always use the "default" queue.
 - `nodes` is the number of machines / computers you request to execute your job. For spike, we always perform execution of a job on one node.
 - `time` is the maximal execution time of your program. By default it is set to 59 minutes and 59 seconds.
 - `mem` is the main memory you request for your job. The scheduler will kill your program if you exceed `time` or `memory` for a significant amount. But keep in mind, too greedy resource definitions will delay execution of your job.
 - `ppn` is the number of compute cored / CPUs for your job. Make sure your program can actually use multiple cores if you request them, otherwise execution is inefficient and other users will start complain about your bad resource estimates!
 
Second entry in the `cluster.json` file is for the rule "demultiplex". This rule will inherit `account`, `queue`, and `nodes` from the "__default__", but override the requested memory to 16GB, the number of cores to just 4 and increase runtime to nearly 2 hours.

## `--cluster "qsub -A {cluster.account} -q {cluster.queue} -l select={cluster.nodes}:ncpus{cluster.ppn}:mem={cluster.mem} -l walltime={cluster.time}"`
You make use of all those above settings by invoking `snakemake` with the flag `--cluster-config cluster.json`, but you also have to define which grid command shall be used to actually submit a job to the grid. For the HPC and the cluster.json file of spike, it looks like `--cluster "qsub -l select={cluster.nodes}:ncpus{cluster.ppn}:mem={cluster.mem} -l walltime={cluster.time}"`. I might recognize the variable names from cluster.json appear here in curly brackets, i.e. those strings will be replaced by the values defined in the cluster.json file.

## `--cluster-status scripts/barnacle_status.py`
From Snakemakes documentation `Status command for cluster execution. This is only considered in combination with the –cluster flag. If provided, Snakemake will use the status command to determine if a job has finished successfully or failed. For this it is necessary that the submit command provided to –cluster returns the cluster job id. Then, the status command will be invoked with the job id. Snakemake expects it to return ‘success’ if the job was successfull, ‘failed’ if the job failed and ‘running’ if the job still runs.` We are using the script https://github.com/sjanssen2/spike/blob/master/scripts/barnacle_status.py for this purpose. (Side note: "barnacle" is the cluster system I was using in San Diego and I was lazy and copy & pasted this script from https://github.com/biocore/oecophylla Thanks Jon for working that one out!)

## `--max-status-checks-per-second 1`
`Snakemake` need to "ping" the scheduler frequently to ask for the status of its jobs. In order to avoid too much asking overhead, I am limiting the number of questions to just one per seconds. This works fine, since `spike` jobs usually run for hours and thus this is no real delay for executing the whole pipeline.

## `--latency-wait 900`
Jobs are executed on specific machines and results will be written in some file. The grid file system than needs to make sure that this file is synchronized with all other machines before execution of a dependent job. This process sometimes takes some time. In my experience, it doesn't hurt to be patient and wait for 900 seconds. If the file does not appear during this long period of time, `snakemake` will treat the job as failed, even though the correct result might pop up later.

## `--use-conda`
HPC admins encourage you to **not** use conda installed packages and instead use their optimized software versions which can be loaded via `module` http://modules.sourceforge.net/
However, since reproducibility is of paramount importance for `spike`, I decided against HPCs suggestion and used many different conda environments to enforce exact program versions that can also easily installed by collaboration partners on their machines. Via the parameter `--use-conda` you enable `snakemake`s fantastic mechanism that automatically downloads and creates conda environments for different rules.

## `-j 100`
This parameter specifies how many of your maybe ten thousands of jobs are submitted at the same time to the scheduler. Don't use much higher numbers as there is the risk that it crashes or slows down the scheduler; not only for you but for *all* users of the HPC!

## `--keep-going`
It might happen that single rules / programs of your pipeline execution fails. Since with `spike` you typically process a multitude of independent samples / trios you don't want to immediately stop execution of all jobs if one fails for one sample. Once you identified and fixed the issue with the one failing job, just re-execute the snakemake command and it will continue from there. Very convenient.

## `-p` and `-r`
The flag `-p` will print the (shell) commands to be executed. Useful for debugging. The flag `-r` reports why a specific rule for a specific input need to be (re)executed. Using both aids debugging and understanding the flow of your pipeline.

# Executing
  1. You should start a new screen https://linuxize.com/post/how-to-use-linux-screen/ session: `screen -S spike`.
  2. Since the login node has very limited resources, you should start an new interactive gird job, with medium memory, one node and one core, but relatively long runtime: `qsub -I -A ngsukdkohi -l mem=10GB -l walltime="40:59:50,nodes=1:ppn=1"`
  3. navigate to your `spike` clone and make sure you configure `spike` correctly (docs to be come), specifically `snupy` credentials and selection of samples to be processed.
  4. Trigger a dry run of the pipeline by using the `-n` flag of snakemake and check that everything looks good: ```snakemake --cluster-config cluster.json --cluster "qsub -A {cluster.account} -q {cluster.queue} -l select={cluster.nodes}:ncpus{cluster.ppn}:mem={cluster.mem} -l walltime={cluster.time}" -j 100 --latency-wait 900 --use-conda --cluster-status scripts/barnacle_status.py --max-status-checks-per-second 1 --keep-going -p -r -n```
  4b. (During development of `spike` it happened to me that `snakemake` wanted to re-execute long running programs because of changes in the Snakefiles, but it would produce identical results. To avoid the waste of compute, I am sometimes "touching" output files to update the time stamp such that snakemake will not reinvoke execution. This harbours the risk of computing with outdated or even incomplete intermediate results! Be careful: replace `-n` with `--touch`.
  5. Trigger actual snakemake run, by removing `-n` (and `--touch`) from the above command.
