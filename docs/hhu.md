## Read this first
`Snakemake` greatly supports execution of a pipeline on a grid compute system with very little overhead. However, it will help reading Snakemake's own documentation to get familiar with with basic concepts and commands: https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution

We are here listing our best practice experience running `spike` on the `High Performance Computing` (HPC) system at University of Duesseldorf: https://www.zim.hhu.de/high-performance-computing.html Great people support users on the HPC via https://rocketchat.hhu.de/ or by phone or email and it is definitively worth reading their Wiki https://wiki.hhu.de/display/HPC/Wissenschaftliches+Hochleistungs-Rechnen+am+ZIM **before** submitting thousands of jobs to their grid.

## `--cluster-config`
Resources on shared computer clusters are not infinit and the job of the admins in to ensure that multiple users get their fair share. This is only possible, if each user behaves nicely, i.e. spends some efforts to find small upper boundaries for their programs needs in terms of memory, CPUs/cores and execution time. Naturally, you are greedy, but admins counteract by putting greedy jobs late in the queue, i.e. you have to wait longer until your greedy job is executed compared to a small job.

`spike` has ~70 different rules, each rule running a different command with different resource demands, depending on typical input sizes. Thus, finding working upper bounds required some experience, which is now distilled in the file https://github.com/sjanssen2/spike/blob/master/cluster.json First entry is the default, which will be applied to every rule if no further matching definition can be found. 
 - In this default rule, we define the `account` to which the used service units get "billed" to. I once had to contact HPC admins and request to create this new account for spike.
 - The HPC has several queues into which your jobs can be submitted. For `spike`, we always use the "default" queue.
 - `nodes` is the number of machines / computers you request to execute your job. For spike, we always perfom execution of a job on one node.
 - `time` is the maximal execution time of your program. By default it is set to 59 minutes and 59 seconds.
 - `mem` is the main memory you request for your job. The scheduler will kill your program if you exceed `time` or `memory` for a significant amount. But keep in mind, too greedy resource definitions will delay execution of your job.
 - `ppn` is the number of compute cored / CPUs for your job. Make sure your program can actually use multiple cores if you request them, otherwise execution is inefficient and other users will start complain about your bad resourse estimates!
 
Seconda entry in the `cluster.json` file is for the rule "demultiplex". This rule will inherit `account`, `queue`, and `nodes` from the "__default__", but override the requested memory to 16GB, the number of cores to just 4 and increase runtime to nearly 2 hours.

You make use of all those settings by invoking `snakemake` with the flag `--cluster-config cluster.json`.
