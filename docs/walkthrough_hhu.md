You can use this document to be guided through the multiple steps to process new (HiSeq) sequencing data.

 0. **Get sample sheet from wetlab crew and validate**
    Unfortunately, people are very lazy when it comes to the sample sheet. Without this document, you don't know what the samples to be sequences are about and cannot proceed processing them. It needs some social engineering to make this point clear. Furthermore, consistent typing of e.g. project names is not always achived. Thus, you need to put extra efforts into sample sheet curation. Additionally, you have to fill in the information about the `spike_entity_id`, `spike_entity_role` and maybe other important data that impact `spike`s processing. [More about sample sheets.](samplesheet.md)
    
    Finally, you should add the new sample sheet (one per flowcell) to the `data/samplesheets` sub-directory (double check it does not contain personal information like patient names!!) and copy it to the input directory of your `spike` installation at the HPC, most likely `/gpfs/project/projects/spike/Inputs/SampleSheets/` (but this can be configures in `config.yaml`).
    
    One erroneously formatted sample sheet can prevent processing of **everything** within `spike`, since *all* sample sheets are read and merged in the very beginning. There is a python function that can validate some aspects - you better use it, see: [More on sample sheets](samplesheet.md) or be very careful with actual execution.
 1. **Transfer raw data to HPC**
    Currently, the HiSeq Instrument is controlled by a Windows PC. The directory `/data/pipeline_in/` of the pipeline server is mounted to this Windows PC as drive `Z:` and new raw data are copied not only to local drives `D:` and `E:`, but also to `Z:`.
    
    SSH into the pipeline server (`134.99.133.111`) and use `rsync` to transfer the raw data from the UKD to the HPC network, best by:
    1. opening a new `screen` session (`screen -S transfer`)
    2. navigate (`cd`) to `/data/pipeline_in/`
    3. `rsync -p -L -r -a -v -e ssh 190327_SN737_0463_BCCN4KACXX jansses@hpc.rz.uni-duesseldorf.de:/gpfs/project/projects/spike/Inputs/Raw_Illumina/`
       
       Of course, you have to replace `190327_SN737_0463_BCCN4KACXX` by the directory of your flowcell, `jansses` by your username at the HPC, and you maybe have to adopt `/gpfs/project/projects/spike/Inputs/Raw_Illumina/` to the actual file path where you want to copy the data to.
       You can also provide multiple directories, in case two flowcells were sequenced in one run.
    4. now wait until transfer is complete. You can disconnect the screen session by hitting the keys `CTRL + a` and than `d` and later reconnect by `screen -dr transfer`.
 2. **Backup raw data to NAS**
    Our backup strategy is quite primitive, as we use multiple individual NAS systems. They are stored in office 13.43.-01.65. The latest one is located in the wet lab room with the HiSeq as "the backup NAS" with IP `10.2.5.12`. We bundle and compress all raw files of each flowcell (which are automatically copied from the controller PC to drive Z:, which is `/data/pipeline_in` on the pipeline server, via this pipeline.
  
    SSH into the pipeline server (`134.99.133.111`)
    1. open a new `screen` session (`screen -S backup`)
    2. navigate (`cd`) to your installation directory of `spike`.
    3. limit the samples known to the pipeline to those of the latest run:
       Information about samples is stored in a [pandas DataFrame](https://pandas.pydata.org/pandas-docs/stable/getting_started/10min.html) called `SAMPLESHEETS` of the main [Snakefile](../Snakefile#L31). You should add a new line above the statement that prints information about the number of samples and projects which starts as `print("%i samples in %i projects."` to subset `SAMPLESHEETS`. For example, you could limit to only those samples of a specific `run` (the better term would have been flow cell, since you can have two flow cells per run - but I haven't had the time to chance that term in spike): `SAMPLESHEETS = SAMPLESHEETS[SAMPLESHEETS['run'].isin(['190327_SN737_0463_BCCN4KACXX'])]`
     4. double check that following settings in file [`config.yaml`](../config.yaml) are properly set:
        1. `dirs`: `prefix` should be the main working data directory for spike. At the pipeline server, this is currently `/data/Spike_data/`
        2. ensure the default temporary directory has sufficient free disk space. This is currently not the case on the pipeline server. Therefore, you must set `dirs`: `tmpdir` to `/data/tmp/`.
        3. check user credentials for the backup NAS. On the pipeline server, it should look like the following, but the password `thisisasecret` must be replaced by the right one.
          ```credentials:
               backup:
                 host: '10.2.5.12'
                 username: "SeqUser"
                 password: "thisisasecret"
                 targetdirectory: "array1/Sequencing_Backups/Illumina_HiSeq"```
      5. execute snakemake, first as a dry run: `snakemake -p --use-conda --cores 30 -r -n`
      6. double check if number of reported samples match your expectations. If so, start the actual processing by `snakemake -p --use-conda --cores 30 -r`, i.e. the same as above without `-n`. Expected runtime for a flow cell backup is roughly one day!
      7. after sucessful execution, `spike` should send a report via email. Once double checked, you can forward this email to the wet lab crew to let them know that the data are savely stored in our backup. They can than free up disk space on the controller PC to prepare the next run. Otherwise, limited hard disk space will not allow to start another run.
      8. From your Windows PC within the UKD network, you should be able to access the web interface of the backup NAS, by entering `https://10.2.5.12` in your favorite browser. User name is `admin`, password is known by e.g. Ute Fischer. Within this interface, you can check free capacity of the NAS. If space is running out, order a new one via a "Bestellschein". You should always have a spare NAS on the shelf in our office, such that you can directly replace. Currently, there are two spare ones - you should be good for ~1 year.
 3. **Subset samples to be processed in Snakefile**
 4. **Execute spike**
    1. send demux- and firstbase report to wetlab crew
 5. Handle errors, monitor progress, re-start failing jobs**
 6. Upload to Snupy
    1. send call report to investigators
