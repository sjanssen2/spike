from collections import OrderedDict


def _format(value):
    try:
        return '%.1f' % int(value)
    except ValueError:
        return value


def annotate(fh_input, fh_output):
    with open(fh_input, 'r') as R:
        with open(fh_output, 'w') as W:
            for line in R.readlines():
                if line.startswith('#'):
                    W.write(line)
                    if line.startswith('##FORMAT=<ID=NR'):
                        W.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Descri'
                                'ption="Number of reads covering variant loca'
                                'tion in this sample",TO_INDEX=0>\n')
                        W.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Descri'
                                'ption="Number of reads covering variant loca'
                                'tion in this sample",FROM_INDEX=0>\n')
                    if line.startswith('##FORMAT=<ID=NV,'):
                        W.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Descri'
                                'ption="Number of reads containing variant in'
                                ' this sample",TO_INDEX=1>\n')
                else:
                    fields = line.rstrip().split('\t')

                    # prepare data in last two columns for easy access
                    data = OrderedDict(
                        [(k, v)
                         for k,v
                         in zip(fields[-2].split(':'), fields[-1].split(':'))])

                    # do the actual computation
                    data.update({'DP': data['NR'].split(',')[0]})
                    parts_nv = data['NV'].split(',')
                    data.update({
                        'AD':
                        ','.join([str(int(data['DP']) - int(parts_nv[0]))] \
                            + parts_nv)})

                    # strang formatting as Sebastian's vcf-mapper does
                    for k in ['GOF']:
                        data.update({k: _format(data[k])})
                    for k in ['GL']:
                        data.update({k: ','.join(map(
                            lambda x: _format(x), data[k].split(',')))})

                    # convert OrderedDict back into VCF format
                    fields[-2] = ':'.join(data.keys())
                    fields[-1] = ':'.join(data.values())

                    W.write('\t'.join(fields)+'\n')
