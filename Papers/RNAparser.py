from Bio import SeqIO

recs = SeqIO.parse(open('/mnt/d/Codes/DataLists/lncRNA/GRCh38_latest_rna.gbff'), 'gb')

fw = open('/mnt/d/Codes/DataLists/lncRNA/len_gc.csv', 'w')

for rec in recs:
    ls = [x for x in rec.__dir__() if not x.startswith('_')]
    name = rec.annotations['molecule_type']
    if name == 'mRNA':
        print(rec.features)
        for feature in rec.features:
            if feature.type == 'CDS':
                break
    else: # RNA, rRNA
        for feature in rec.features:
            if feature.type == 'ncRNA':
                break
    start = feature.location.start
    end = feature.location.end
    seq = rec.seq[start: end]
    gc = (seq.count('G') + seq.count('C')) / len(seq)
    length = abs(end - start)
    fw.write('%s,%d,%.5f\n' % (name, length, gc))
    break
fw.close()
