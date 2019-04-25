The samplesheet is the place where all information about samples and sequencing runs are accumulated - thus, it is the most valuable piece of data for our projects. Sequence data are costly to generate, but worthless without knowing what samples have been processed. So - please - put some love into curation of your samplesheets!

## Structure
The samplesheets used in `spike` are based on the file format of the demultiplexing sheets used with [Illumina's bcl2fastq software](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq_letterbooklet_15038058brpmi.pdf), i.e. comma separated lines.
The file consists of 4 sections: `Header`, `Reads`, `Settings`, `Data`.
#### Header
The header section has multiple lines consisting of two columns. The first column is the *key*, the second the *value*. (There can be more empty columns indicated by multiple `,` characters)  Ordering of line is arbitrary.
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
#### Reads
This small section reflects the `kind_of_run`. I am not confident if it impacts `bcl2fastq`'s behaviour, but you better keep it in sync with `kind_of_run`, i.e. one line for single end, two lines for paired end holding the number of sequenced cycles, most likely `101`.
#### Settings
I honestly don't understand the meaning of it (please update my knowledge gap here) but I see that it is always identically. Thus, I suggest you copy and paste as is:
```
[Settings]
ReverseComplement,0
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
#### Data
This is the most important section as it is used to demultiplex the raw data and controls the behaviour of the `spike` pipeline.
The inherited field from Illumina's demultiplexing sheets are:
 - `Lane`: a HiSeq flowcell has 8 separate lanes and demultiplexing only occures within one lane. MiSeq machines only have 1 lane. Leave the lane number empty when dealing with situations where you don't have raw data, but only fastq files. This is the case when sequencing was performed by Macrogen.
 - `Sample_ID`: this is the name of the sample. Note that sequencing data is merged by `spike` if identical `Sample_ID`s (`Sample_Project` must also be identical) are found either in different lanes (that is a typical setup for `Keimbahn` samples) or on different flowcells (if a sample was re-sequenced to obtain a better coverage). 
 - `Sample_Name`: needs to be empty!
 - `Sample_Plate`: is ignored by `bcl2fastq` and `spike`
 - `Sample_Well`: is ignored by `bcl2fastq` and `spike`
 - `I7_Index_ID`: is the identifyer of the barcode sequence
 - `index`: is the barcode sequence used to demultiplex. Ensure you don't have any typos in here as this would break the sample demultiplexing!
 - `Sample_Project`: Samples are organized into different, independent projects. Processing, data delivery and access rights within Snupy happens project-wise. Thus, it is important that you type in existing, correct names.
 - `Description`: feel free to leave notes, this field is ignored otherwise.

Additional field for `spike`:
The following field fundamentally control `spike`s type of processing. Thus, double check everything is correct here.
 - `spike_entity_id`: The "name" of the family in case of the Keimbahn project, e.g. `KB0011`. Or e.g. `mTAS_360` for the mouse project. This groups a set of samples for different use cases, either being the plain demultiplexing if the field remains empty, trio computation for father-mother-child samples or tumornormal substraction for tumor-healthy tissue samples. For more see below:
 - `spike_entity_role`: Within a `spike_entity_id` a sample has a specific "role". For the "trio" use case, there need to be at least the three roles `father`, `mother`, and `patient` (which is the child). An additional role is `sibling` for whome another trio computation is performed as for the patient. The second use case is "tumornormal" where a `spike_entity_id` needs to define `tumor` and `healthy` (= normal) sample. Edge cases are roles like `tumor_father` or `tumor_patient` where the same sample is not only used for "tumornormal" computation, but also for trio computation.
 - `spike_ignore_sample`: you might want to completely exclude a sample for further processing if you e.g. know that sample preparation failed for biological reasons. Thus, you don't have to run into failing processing steps again and again.

Optinal fields for `spike`:
 - `spike_notes`: there are some badly designed combinations of barcodes (mixing amplicons and WES) where you need to lower `bcl2fastq`'s allowed error number to succeede with demultiplexing (with the risk of barcode hopping). You can mark all samples of such a lane by putting `colliding barcode` in this field.
 - `capture_kit`: typically, the `Sample_Project` defines the WES capture kit, but there are a few instances where the kit is changed within the same project. Thus, you can define the kit here sample-wise instead of project-wise.
 
## Global samplesheets table
In the beginning of any `Spike` activity, first the configuration file is read, than **all** samplesheets in the defined directory are read and combined. Thus, violating the file format for one samplesheet can break everything! Thus, I suggest you copy and paste an existing samplesheet - edit the information as needed and save as a new file. Only if you really really know what you are doing you can create a samplesheet from scratch. Working examples are located in the directory `data/samplesheets/`.

Furthermore, there is a python function [validate_samplesheet](https://github.com/sjanssen2/spike/blob/ea474be625d9d8550dcf4be3950a25806120c4cd/scripts/parse_samplesheet.py#L86) which you should execute on a new samplesheet to programmatically check for common mistakes, like typos in Sample_IDs spread across multiple lanes, use of restricted character like `-`, notification about missing columns, warning about empty `Sample_Project`, untypical `spike_entity_roles` or mismatching `Sample_ID` to `spike_entity_id`, ... Just fire up a [jupyter notebook](https://jupyter.org/) and run the validation function!
