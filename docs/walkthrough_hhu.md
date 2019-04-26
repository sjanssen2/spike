You can use this document to be guided through the multiple steps to process new (HiSeq) sequencing data.

 0. **Get sample sheet from wetlab crew and validate**
    Unfortunately, people are very lazy when it comes to the sample sheet. Without this document, you don't know what the samples to be sequences are about and cannot proceed processing them. It needs some social engineering to make this point clear. Furthermore, consistent typing of e.g. project names is not always achived. Thus, you need to put extra efforts into sample sheet curation. Additionally, you have to fill in the information about the `spike_entity_id`, `spike_entity_role` and maybe other important data that impact `spike`s processing. [More on sample sheets](samplesheet.md)
    
    Finally, you should add the new sample sheet (one per flowcell) to the `data/samplesheets` sub-directory (double check it does not contain personal information like patient names!!) and copy it to the input directory of your `spike` installation at the HPC, most likely `/gpfs/project/projects/spike/Inputs/SampleSheets/` (but this can be configures in `config.yaml`).
    
    One erroneously formatted sample sheet can prevent processing of anything within `spike`, since *all* sample sheets are read and merged in the very beginning. There is a python function that can validate some aspects - you better use it, see: [More on sample sheets](samplesheet.md) or be very careful with actual execution.
 1. **Transfer raw data to HPC**
    Currently, the HiSeq Instrument is controlled by a Windows PC. The directory `/data/pipeline_in/` of the pipeline server is mounted to this Windows PC as drive `Z:` and new raw data are copied not only to local drives `D:` and `E:`, but also to `Z:`.
    
    SSH into the pipeline server (`134.99.133.111`) and use `rsync` to transfer the raw data from the UKD to the HPC network, best by:
    1. opening a new `screen` session (`screen -S transfer`)
    2. navigate (`cd`) to `/data/pipeline_in/`
    3. `rsync -p -L -r -a -v -e ssh 190327_SN737_0463_BCCN4KACXX jansses@hpc.rz.uni-duesseldorf.de:/gpfs/project/projects/spike/Inputs/Raw_Illumina/`
       
       Of course, you have to replace `190327_SN737_0463_BCCN4KACXX` by the directory of your flowcell, `jansses` by your username on the HPC, and you maybe have to adopt `/gpfs/project/projects/spike/Inputs/Raw_Illumina/` to the actual file path where you want to copy the data to.
       You can also provide multiple directories, in case two flowcells were sequenced in one run.
    4. now wait until transfer is complete. You can disconnect the screen session by hitting the keys `CTRL + a` and than `d` and later reconnect by `screen -dr transfer`.
 2. **Backup raw data to NAS**
    *1. order / configure new NAS if current one is full*
    *2. let wetlab crew know that they can delete raw data from sequencing PC.*
 3. **Subset samples to be processed in Snakefile**
 4. **Execute spike**
    1. send demux- and firstbase report to wetlab crew
 5. Handle errors, monitor progress, re-start failing jobs**
 6. Upload to Snupy
    1. send call report to investigators
