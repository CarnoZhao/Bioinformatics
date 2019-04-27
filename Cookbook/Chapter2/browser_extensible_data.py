from collections import defaultdict
import re
import HTSeq

class Browser_Entensible_Data():
    '''
    using HTseq python package to deal with bed file
    '''

    def __init__(self, filename = '../../../DataLists/CookbookData/LCT.bed'):
        self.bed = HTSeq.BED_Reader(filename)

    def print_information(self):
        feature_types = defaultdict(int)
        for rec in self.bed:
            last_rec = rec
            feature_types[re.search('([A-Z]+)', rec.name).group(0)] += 1
        print('feature_types: ', feature_types)
        print('record: ', last_rec)
        print('record name: ', last_rec.name)
        print('type of record: ', type(last_rec))
        interval = last_rec.iv
        print('record interval: ', interval)
        print('type of interval: ', type(interval))
        print('chromosome of interval: ', interval.chrom)
        print('start and end of interval: ', interval.start, interval.end)
        print('interval strand: ', interval.strand)
        print('length of interval: ', interval.length)
        print('interval.start_d: ', interval.start_d, ' regarless of the strand, can be seen as the absolute position of the gene in a chromosome')
        print('interval.start_as_pos and its type', interval.start_as_pos, type(interval.start_as_pos))

    def CDS_statistics(self):
        '''

        '''
        exon_start = None
        exon_end = None
        sizes = []
        for rec in self.bed:
            if not rec.name.startswith('CCDS'):
                continue
            interval = rec.iv
            exon_start = min(interval.start, exon_start or interval.start)
            exon_end = max(interval.end, exon_end or interval.end)
            sizes.append(interval.length)
        sizes.sort()
        print('Num exons: %d / Begin: %d / End %d' % (len(sizes), exon_start, exon_end))
        print('Smaller exon: %d / Larger exon: %d / Mean size %.1f' % (sizes[0], sizes[-1], sum(sizes) / len(sizes)))

if __name__ == '__main__':
    bed = Browser_Entensible_Data()
    # bed.print_information()
    bed.CDS_statistics()

