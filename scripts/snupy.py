import hashlib
import pandas as pd
import requests
from requests.auth import HTTPBasicAuth
import gzip
import shutil
import tempfile
from os.path import join
import sys
import time
import re
from scripts.parse_samplesheet import get_role


def get_md5sum(fname):
    hash_md5 = hashlib.md5()
    with gzip.open(fname, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def get_toolname_from_stepname(config, filename):
    if config['stepnames']['platypus_filtered'] in filename:
        return 'Platypus'
    elif config['stepnames']['gatk_CombineVariants'] in filename:
        return 'GATK'
    elif config['stepnames']['mutect'] in filename:
        return 'Mutect'
    elif config['stepnames']['merge_somatic'] in filename:
        return 'VarScan2'
    elif config['stepnames']['writing_headers'] in filename:
        return 'VarScan2'
    elif config['stepnames']['excavator_somatic'] in filename:
        return 'Excavator2'
    elif config['stepnames']['excavator_trio'] in filename:
        return 'Excavator2'
    else:
        raise ValueError("Unexpected tool.")


def check_snupy_status(request_result):
    status = request_result.headers.get('status')
    if (status is None) and \
       ('<title>401 Unauthorized</title>' in request_result.text):
        raise ValueError(("Snupy password is not correct. Please check what "
                          "you have provided in the config.yaml file."))
    assert(status == '200 OK')


def is_alias_sample(filename, samplesheets):
    # check if this sample is used as an alias ...
    fastq_prefix = '/'.join(filename.split('/')[-2:]).split('.')[0]
    return samplesheets[(samplesheets['fastq-prefix'] == fastq_prefix) & (samplesheets['is_alias'] == True)]


def get_snupy_sample_name(project, entity, filename, config, samplesheets, _type):
    name = ''

    # check if this sample is used as an alias ...
    aliases = is_alias_sample(filename, samplesheets)
    if aliases.shape[0] > 0:
        # ... yes, we have a sample alias at hand
        filename = filename.replace(aliases['fastq-prefix'].iloc[0], '%s/%s' % (project, aliases['Sample_ID'].iloc[0]))

    # original sample name if background else spike_entity_id
    sample_name = entity
    if _type == 'background':
        sample_name = filename.split('/')[-1].split('.')[0]

    # IF sample has been sequence on more than one flowcell (i.e. at different rundates)
    # the dates are added to the name to give the user a hint about this fact
    runs = samplesheets[(samplesheets['Sample_Project'] == project) &
                        ((samplesheets['spike_entity_id'] == entity) | (samplesheets['Sample_ID'] == entity)) &
                        ((samplesheets['Sample_ID'] == sample_name) | (_type != 'background'))]['run'].unique()
    name += '+'.join(sorted(runs))
    name += '_' + project
    name += '/' + sample_name

    # snv type
    snvtype = None
    if get_toolname_from_stepname(config, filename) == 'GATK':
        snvtype = 'snp_indel'
    elif get_toolname_from_stepname(config, filename) == 'Platypus':
        snvtype = 'indel'
    elif get_toolname_from_stepname(config, filename) == 'VarScan2':
        if config['stepnames']['merge_somatic'] in filename:
            snvtype = 'somatic'
        elif config['stepnames']['writing_headers'] in filename:
            snvtype = 'denovo'
    elif get_toolname_from_stepname(config, filename) == 'Mutect':
        snvtype = 'somatic'
    elif get_toolname_from_stepname(config, filename) == 'Excavator2':
        if config['stepnames']['excavator_trio'] in filename:
            snvtype = 'trio.cnv'
        elif config['stepnames']['excavator_somatic'] in filename:
            snvtype = 'somatic.cnv'
        else:
            raise ValueError("Unknown SNV type for Excavator2")
    else:
        raise ValueError("Unknown SNV type")
    name += ".%s" % snvtype

    # program
    name += ".%s" % get_toolname_from_stepname(config, filename).lower().replace('platypus', 'ptp').replace('varscan2', 'varscan')

    # file ending
    #name += '.vcf'

    return name


def get_snupy_parser(config, filename):
    if config['stepnames']['excavator_somatic'] in filename:
        return 'VcfFileExcavator'
    if config['stepnames']['excavator_trio'] in filename:
        return 'VcfFileExcavator'
    if (config['stepnames']['merge_somatic'] in filename) or (config['stepnames']['writing_headers'] in filename):
        return 'VcfFileVarscan'
    return 'VcfFile'


def get_upload_content(project, entity, input, config, samplesheets, tmpdir, _type, snupy_instance):
    data = pd.DataFrame(index=input)
    zippedfiles = []
    for i, file in enumerate(data.index):
        sys.stderr.write('compressing %i of %i files: %s\n' % (i+1, len(data.index), file.split('/')[-1]))
        with open(file, 'rb') as f_in:
            file_gz = join(tmpdir, file.split('/')[-1])+'.gz'

            aliases = is_alias_sample(file, samplesheets)
            with gzip.open(file_gz, 'wb') as f_out:
                #shutil.copyfileobj(f_in, f_out)
                for line in f_in.readlines():
                    # those are the first mandatory 8 columns:
                    # https://en.wikipedia.org/wiki/Variant_Call_Format
                    if (aliases.shape[0] > 0) and (line.startswith(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")):
                        # search for the original sample name and replace by alias
                        line = line.replace(bytes(aliases['fastq-prefix'].iloc[0].split('/')[-1], encoding='utf-8'),
                                            bytes(aliases['Sample_ID'].iloc[0], encoding='utf-8'))
                        sys.stderr.write("  rename '%s' --> '%s'\n" % (aliases['fastq-prefix'].iloc[0].split('/')[-1], aliases['Sample_ID'].iloc[0]))
                    f_out.write(line)
            zippedfiles.append(file_gz)
    data['zipped'] = zippedfiles
    data['contact'] = config['credentials']['snupy'][snupy_instance]['username']
    data['institution_id'] = config['snupy_ids'][snupy_instance]['institues']['UKD']
    data['tags[TOOL]'] = list(map(lambda x: config['snupy_ids'][snupy_instance]['tools'][get_toolname_from_stepname(config, x)], data.index))
    data['type'] = list(map(lambda x: get_snupy_parser(config, x), data.index))
    data['organism_id'] = config['snupy_ids'][snupy_instance]['organisms'][config['projects'][project]['species']]
    data['md5checksum'] = list(map(get_md5sum, data['zipped'].values))
    data['name'] = list(map(lambda x: get_snupy_sample_name(project, entity, x, config, samplesheets, _type), data.index))

    payload = dict()
    files = dict()
    for i, (filepath, row) in enumerate(data.iterrows()):
        for col in row.index:
            if col == 'zipped':
                files['vcf_file_templates[%i][content]' % (i+1)] = (row['zipped'].split('/')[-1], open(row['zipped'], 'rb'), 'application/gzip', {'Expires': '0'})
            elif col == 'role':
                continue
            else:
                payload['vcf_file_templates[%i][%s]' % (i+1, col)] = row[col]

    return payload, files, data


def upload_to_snupy(project, entity, input, config, samplesheets, output, log, _type, snupy_instance):
    tmpdir = tempfile.mkdtemp()

    payload, files, data = get_upload_content(project, entity, input, config, samplesheets, tmpdir, _type, snupy_instance)
    # print(payload)
    # print(files)
    # return None
    r = requests.post(
        str(config['credentials']['snupy'][snupy_instance]['host'] + '/vcf_files/batch_submit.json'),
        data=payload,
        auth=HTTPBasicAuth(config['credentials']['snupy'][snupy_instance]['username'], config['credentials']['snupy'][snupy_instance]['password']),
        verify=False,
        files=files
        )

    with open(log[0], 'w') as f:
        f.write('==== headers ====\n')
        f.write(str(r.headers))
        f.write("\n\n")
        f.write('==== text ====\n')
        f.write(str(r.text))

    # assert HTTP communication successful
    check_snupy_status(r)
    shutil.rmtree(tmpdir, ignore_errors=True)

    # assert VCF upload to Snupy was successful
    if len(r.json()['error']) != 0:
        raise ValueError('Snupy upload failed with error(s): "%s"\n' % '" and "'.join({err['output']['notice'] for err in r.json()['error']}))
    # parse response from snupy into pandas table
    res = pd.DataFrame(r.json()['success'])[["id", "sample_names", "status", "name", "updated_at", "created_at", "md5checksum"]]
    res.columns = map(lambda x: 'snupy_%s' % x, res.columns)

    # add data from spike
    res['spike_entity_id'] = entity
    res['spike_project'] = project
    res['spike_type'] = _type

    res = res.merge(data[['md5checksum', 'tags[TOOL]']], left_on='snupy_md5checksum', right_on='md5checksum')

    # write snupys result table as output file
    res.to_csv(str(output[0]), sep="\t")

    return res


def extractsamples(uploadtable, config, samplesheets, output, log, _type, snupy_instance):
    uploaded = pd.read_csv(uploadtable[0], sep="\t", index_col=0, dtype=str)
    uploaded.rename(columns={'snupy_id': 'vcf_file_id'}, inplace=True)

    extracted = uploaded.copy()
    extracted['tags'] = None
    extracted['min_read_depth'] = 5
    extracted['specimen_probe_id'] = " "
    extracted['contact'] = config['credentials']['snupy'][snupy_instance]['username']
    extracted['ignorefilter'] = '2'
    extracted['filters'] = 'PASS'
    extracted['info_matches'] = ""
    cache_users = dict()

    for idx, row in extracted.iterrows():
        if _type == 'background':
            extracted.loc[idx, 'tags'] = '{"DATA_TYPE":"snv and indel"}'
            extracted.loc[idx, 'snupy_Samples'] = row['snupy_name'].split('/')[-1].split('.')[0]
            if samplesheets[(samplesheets['Sample_ID'] == extracted.loc[idx, 'snupy_Samples']) &
                            (samplesheets['Sample_Project'] == row['spike_project']) &
                            (samplesheets['spike_entity_id'] == row['spike_entity_id'])]['spike_entity_role'].iloc[0] in ['healthy', 'father', 'mother', 'sibling']:
                extracted.loc[idx, 'ignorefilter'] = '1'
        elif _type == 'tumornormal':
            extracted.loc[idx, 'tags'] = '{"DATA_TYPE":"somatic"}'
            if row['tags[TOOL]'] == str(config['snupy_ids'][snupy_instance]['tools']['VarScan2']):
                extracted.loc[idx, 'snupy_Samples'] = 'TUMOR'
            elif row['tags[TOOL]'] == str(config['snupy_ids'][snupy_instance]['tools']['Mutect']):
                extracted.loc[idx, 'snupy_Samples'] = get_role(row['spike_project'], row['spike_entity_id'], 'tumor', samplesheets).split('/')[-1]
            elif row['tags[TOOL]'] == str(config['snupy_ids'][snupy_instance]['tools']['Excavator2']):
                extracted.loc[idx, 'snupy_Samples'] = get_role(row['spike_project'], row['spike_entity_id'], 'tumor', samplesheets).split('/')[-1]
                extracted.loc[idx, 'min_read_depth'] = 0
                extracted.loc[idx, 'tags'] = '{"DATA_TYPE":"cnv"}'
        elif _type == 'trio':
            extracted.loc[idx, 'tags'] = '{"DATA_TYPE":"denovo"}'
            if row['tags[TOOL]'] == str(config['snupy_ids'][snupy_instance]['tools']['VarScan2']):
                extracted.loc[idx, 'snupy_Samples'] = 'Child'
            elif row['tags[TOOL]'] == str(config['snupy_ids'][snupy_instance]['tools']['Excavator2']):
                extracted.loc[idx, 'snupy_Samples'] = get_role(row['spike_project'], row['spike_entity_id'], 'patient', samplesheets).split('/')[-1]
                extracted.loc[idx, 'min_read_depth'] = 0
                extracted.loc[idx, 'tags'] = '{"DATA_TYPE":"cnv"}'
        extracted.loc[idx, 'nickname'] = row['snupy_name'].split('/')[-1]
        extracted.loc[idx, 'project'] = str(config['projects'][row['spike_project']]['snupy'][snupy_instance]['project_id'])
        if row['spike_project'] not in cache_users:
            r = requests.get(config['credentials']['snupy'][snupy_instance]['host'] + ('/experiments/%s.json' % extracted.loc[idx, 'project']),
                auth=HTTPBasicAuth(config['credentials']['snupy'][snupy_instance]['username'], config['credentials']['snupy'][snupy_instance]['password']),
                verify=False)
            cache_users[row['spike_project']] = ';'.join(sorted(map(lambda x: x['name'], r.json()['users'])))
        extracted.loc[idx, 'users'] = cache_users[row['spike_project']]

    extracted.rename(columns={'snupy_vcfID': 'vcf_file_id',
                              'snupy_Samples': 'vcf_sample_name',
                              'snupy_name': 'name',
                              'spike_entity_id': 'patient',
                              }, inplace=True)
    res = extracted[['vcf_file_id', 'vcf_sample_name', 'tags', 'name', 'nickname', 'patient', 'min_read_depth',
                     'project', 'specimen_probe_id', 'contact', 'users', 'ignorefilter', 'filters', 'info_matches']]

    # write sample sheet into temporary file
    _, fp = tempfile.mkstemp()
    res.to_csv(fp, sep="\t", index=False, quotechar="'")

    # uploading the generated sample extraction sheet
    files = {'sample[sample_sheet]': (
        'SampleExtractionSheet.csv',
        open(fp, 'rb'), 'application/txt', {'Expires': '0'})}

    r = requests.post('%s/samples.json' % config['credentials']['snupy'][snupy_instance]['host'],
                      auth=HTTPBasicAuth(config['credentials']['snupy'][snupy_instance]['username'], config['credentials']['snupy'][snupy_instance]['password']),
                      verify=False,
                      files=files,
                     )
    with open(log[0], 'w') as f:
        f.write('==== headers ====\n')
        f.write(str(r.headers))
        f.write("\n\n")
        f.write('==== text ====\n')
        f.write(str(r.text))

    check_snupy_status(r)
    shutil.rmtree(fp, ignore_errors=True)

    # parse snupy's response
    snupy_table = pd.DataFrame(r.json())

    # write snupys result table as output file
    snupy_table.to_csv(str(output[0]), sep="\t")

    if any(map(lambda x: 'error' in x, r.json())):
        sys.stderr.write(snupy_table.to_string())
        raise ValueError('Something went wrong when samples should be extracted in snupy: "%s"' % '" and "'.join({x['error'] for x in r.json()}))
