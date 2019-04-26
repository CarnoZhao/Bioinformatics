from Bio import SeqIO
import gzip

class Fastq_Analysis():
    '''
    some simple analysis of fastq file, including the qualities and the sequence itself
    '''

    def __init__(self, filename = '../../../DataLists/CookbookData/SRR003265.filt.fastq.gz'):
        self.filename = filename
        self.recs = SeqIO.parse(gzip.open(filename, 'rt', encoding = 'utf-8'), 'fastq')
    
    def print_simple_information(self):
        rec = next(self.recs)
        print('rec id:', rec.id, rec.description, rec.seq)
        print(rec.letter_annotations)


if __name__ == '__main__':
    fq = Fastq_Analysis()
    fq.print_simple_information()
    
