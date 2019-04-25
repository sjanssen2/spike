The samplesheet is the place where all information about samples and sequencing runs are accumulated - thus, it is the most valuable piece of data for our projects. Sequence data are costly to generate, but worthless without knowing what samples have been processed. So - please - put some love into curation of your samplesheets!

## Structure
The samplesheets used in `spike` are based on the file format of the demultiplexing sheets used with [Illumina's bcl2fastq software](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf), i.e. comma separated lines.
The file consists of 4 sections: `Header`, `Reads`, `Settings`, `Data`.
#### Header
The header section has multiple lines consisting of two columns. The first column is the *key*, the second the *value*. (There can be more empty columns indicated by multiple `,` characters)
You can process your data without those information, but I strongly suggest to also type in these fields as they ease organization of samples.
Important keys are:
 - `Backup_NAS_name`: We use multiple NAS devices to back up the raw Illumina data (see XXX). To avoid booting all NAS (we currently have 9 and this number is counting) when searching for a backup, it is useful to note down on which of those NAS this flow cell data have been stored. Do youself the favor and enter this tiny piece of information instead of spending hours to find the data if necessary!
 - `kind_of_run`: e.g. `2x101bp` means that the sequencer ran in paired end mode (`2x`) for `101bp` cycles. This is our current default for WES sequencing, but we also have data for single end or longer reads.
 - `Operator`: is the person responsible for prep'ing the DNA, loading the samples onto the flowcell and operating the sequencer. You want to know this to be able to later ask the right person if something looks wired. Currently, these people are "Katayoun Alemazkour", "Daniel Scholtyssik" and "Frauke Meyer" (but this information might be outdated)
 - `Investigator`: is the person waiting for the data, e.g. "Ute Fischer" or "Daniel Picard". You want to notify them once you finished processing. Thus, it is handy to know who your "customer" was.
 - `Date`: the date of starting the run on the Illumina machine, unfortunately in american format, i.e. mm/dd/yy. This should be in sync with the directory name of the raw Illumina data.
 - `Instrument Type`: for now, this is our HiSeq machine, i.e. `HiSeq 15000/2500` or `MiSeq` for 16s data from the BMFZ or something else e.g. if sequencing was outsourced to another provider.
 - `Assay`: the capture kit, like `Agilent_Sure_SelectXT`
 - `Index Adapters`: `Agilent_Sure_SelectXT`
 - `Chemistry`: the chemisty version and lot used for the run, e.g. `20324563 Truseq PE cluster kit V3`
 
