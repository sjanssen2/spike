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
from bs4 import BeautifulSoup
import re


MAP_INSTITUTES = {'UKD': 4}
MAP_ORGANISMS = {'homo sapiens': 1,
                 'mus musculus': 2}
MAP_TOOLS = {'GATK': 118,
             'GATK_RELAXED': 119,
             'GATK_TRIO': 120,
             'Mutect': 121,
             'VarScan2': 122,
             'Platypus': 123,
             'Excavator2': 124,}


def get_md5sum(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
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
    else:
        raise ValueError("Unexpected tool.")


def get_snupy_sample_name(project, entity, filename, config, samplesheets, _type):
    name = 'spike_'

    # original sample name if background else spike_entity_id
    sample_name = entity
    if _type == 'background':
        sample_name = filename.split('/')[-1].split('.')[0]
    name += sample_name

    # IF sample has been sequence on more than one flowcell (i.e. at different rundates)
    # the dates are added to the name to give the user a hint about this fact
    seqdates = samplesheets[(samplesheets['Sample_Project'] == project) &
                            (samplesheets['spike_entity_id'] == entity) &
                            ((samplesheets['Sample_ID'] == sample_name) | (_type != 'background'))]['run'].unique()
    seqdates = list(map(lambda x: x.split('_')[0], seqdates))
    if len(seqdates) > 1:
        seqdates = '_merged%s' % ('+'.join(seqdates))
    else:
        seqdates = ""
    name += seqdates

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
    else:
        raise ValueError("Unknown SNV type")
    name += ".%s" % snvtype

    # program
    name += ".%s" % get_toolname_from_stepname(config, filename).lower().replace('platypus', 'ptp').replace('varscan2', 'varscan')

    # file ending
    #name += '.vcf'

    return name

def get_upload_content(project, entity, input, config, samplesheets, tmpdir, _type):
    data = pd.DataFrame(index=input)
    zippedfiles = []
    for i, file in enumerate(data.index):
        sys.stderr.write('compressing %i of %i files: %s\n' % (i+1, len(data.index), file.split('/')[-1]))
        with open(file, 'rb') as f_in:
            file_gz = join(tmpdir, file.split('/')[-1])+'.gz'
            with gzip.open(file_gz, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            zippedfiles.append(file_gz)
    data['zipped'] = zippedfiles
    data['contact'] = config['credentials']['snupy']['username']
    data['institution_id'] = MAP_INSTITUTES['UKD']
    data['tags[TOOL]'] = list(map(lambda x: MAP_TOOLS[get_toolname_from_stepname(config, x)], data.index))
    data['type'] = list(map(lambda x: 'VcfFileVarscan' if get_toolname_from_stepname(config, x) == 'VarScan2' else 'VcfFile', data.index))
    data['organism_id'] = MAP_ORGANISMS[config['projects'][project]['species']]
    data['md5checksum'] = list(map(get_md5sum, data.index))
    data['name'] = list(map(lambda x: get_snupy_sample_name(project, entity, x, config, samplesheets, _type), data.index))

    payload = dict()
    files = dict()
    for i, (filepath, row) in enumerate(data.iterrows()):
        for col in row.index:
            if col == 'zipped':
                files['vcf_file_templates[%i][content]' % (i+1)] = (row['zipped'].split('/')[-1], open(row['zipped'], 'rb'), 'application/gzip', {'Expires': '0'})
            else:
                payload['vcf_file_templates[%i][%s]' % (i+1, col)] = row[col]


    return payload, files

def upload_to_snupy(project, entity, input, config, samplesheets, output, log, _type):
    tmpdir = tempfile.mkdtemp()

    payload, files = get_upload_content(project, entity, input, config, samplesheets, tmpdir, _type)
    print(payload)
    print(files)
    # return None
    r = requests.post(
        str(config['credentials']['snupy']['host'] + '/vcf_files/batch_submit'),
        data=payload,
        auth=HTTPBasicAuth(config['credentials']['snupy']['username'], config['credentials']['snupy']['password']),
        verify=False,
        files=files
        )

    with open(log[0], 'w') as f:
        f.write('==== headers ====\n')
        f.write(str(r.headers))
        f.write("\n\n")
        f.write('==== text ====\n')
        f.write(str(r.text))

    assert(r.headers.get('status') == '200 OK')
    shutil.rmtree(tmpdir, ignore_errors=True)

    # parse response from snupy into pandas table
    soup = BeautifulSoup(r.text, 'html.parser')
    res_table = soup.find("div", re.compile('fail|success')).table
    cols = [name.text for name in res_table.thead.find_all('th')]
    values = []
    for row in res_table.find_all('tr'):
        values.append([value.text for value in row.find_all('td')])
    res = pd.DataFrame(values, columns=cols).dropna()

    if ('Status' not in res.columns) or (list(res['Status'].unique()) != ['CREATED']):
        sys.stderr.write(res.to_string())
        raise ValueError("Something went wrong when files should be uploaded to snupy.")

    # write snupys result table as output file
    res.to_csv(str(output[0]), sep="\t", index_label=True)

    return res
