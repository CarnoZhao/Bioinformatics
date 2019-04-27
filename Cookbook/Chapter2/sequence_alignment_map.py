import pysam
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

class Sequence_Alignment_Map():
    '''
    Some basic analysis of sam/bam file using pysam package
    '''

    def __init__(self, filename = '../../../DataLists/CookbookData/NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam'):
        self.bam = pysam.AlignmentFile(filename, 'rb')
    
    def bam_header_information(self):
        '''
        print the information of headers of bam file
        '''
        headers = self.bam.header
        for recordtype, records in headers.items():
            print(recordtype)
            for i, record in enumerate(records):
                if i > 4:
                    print('\t......')
                    break
                if type(record) == dict:
                    print('\t%d' % (i + 1))
                    for field, value in record.items():
                        print('\t\t%s\t%s' % (field, value))
                else:
                    print('\t%s' % record)
    
    def bam_record_information(self):
        '''
        print the information of a record (a read) of bam file
        '''
        global rec
        for rec in self.bam:
            if rec.cigarstring.find('M') > -1 and rec.cigarstring.find('S') > -1 and not rec.is_unmapped and not rec.mate_is_unmapped:
                break
        self.simplify_print('rec.query_', ('name',))
        self.simplify_print('rec.reference_', ('id', 'start', 'end'))
        self.simplify_print('', ('self.bam.getrname(rec.reference_id)',))
        self.simplify_print('', ('rec.cigarstring',))
        self.simplify_print('rec.query_alignment_', ('start', 'end', 'length'))
        self.simplify_print('rec.next_reference_', ('id', 'start'))
        self.simplify_print('', ('rec.template_length',))
        self.simplify_print('rec.is_', ('paired', 'proper_pair', 'unmapped'))
        self.simplify_print('', ('rec.mapping_quality',))
        self.simplify_print('rec.query_', ('qualities', 'alignment_qualities', 'sequence'))


    def mapped_position_distribution(self, chrom = '20', start = 0, end = 10000000):
        '''
        plot the distribution of successfully mapped number according to the position
        '''
        counts = [0 for _ in range(76)]
        for n, rec in enumerate(self.bam.fetch(chrom, start, end)):
            for i in range(rec.query_alignment_start, rec.query_alignment_end):
                counts[i] += 1
        fig, ax = plt.subplots(figsize = (16, 9))
        ax.plot(range(1, 77), [x / (n + 1) for x in counts])
        plt.savefig('mapped_position_distribution.pdf')


    def simplify_print(self, prefix, suffix):
        for suf in suffix:
            print(prefix + suf + ': ', eval(prefix + suf), end = '\n')
        print()

if __name__ == '__main__':
    sam = Sequence_Alignment_Map()
    # sam.bam_header_information()
    # sam.bam_record_information()
    sam.mapped_position_distribution()
