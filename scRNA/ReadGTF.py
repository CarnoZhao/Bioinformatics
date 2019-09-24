import pandas as pd

names = []
types = []
chroms = []
with open('D:/Download/Mus_musculus.GRCm38.97.chr.gtf') as gtf:
    for l in gtf:
        if l.startswith('#'):
            continue
        chrom = l.split('\t')[0]
        # if chrom != 'MT':
        #     continue
        # annos = l.split('\t')[-1]
        # annos = l.split('\t')[-1].strip()[:-1].replace('; ', ',"').replace(' ', '":')
        # d = eval('{"' + annos + '}')
        # names.append(d['gene_name'])
        # types.append(d['gene_biotype'])
#         if 'rRNA' not in l:
#             continue
#         else:
#             annos = l.split('\t')[-1]
#             annos = l.split('\t')[-1].strip()[:-1].replace('; ', ',"').replace(' ', '":')
#             d = eval('{"' + annos + '}')
#             names.append(d['gene_name'])
#             types.append(d['gene_biotype'])
# with open('D:/Codes/DataLists/Output/rRNA.txt', 'w') as rRNA:
#     rRNA.write('\n'.join(names))
print('\t'.join(names))