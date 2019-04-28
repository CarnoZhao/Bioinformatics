from Bio import SeqIO, SeqUtils
from collections import defaultdict

class General_Genome():
    '''
    general information about a genome fasta file
    '''
    def __init__(self, filename = '../../../DataLists/CookbookData/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta'):
        self.fasta = SeqIO.parse(filename, 'fasta')
        
    def print_all_description(self):
        '''
        print the description of all fasta records, which is separated to different chromosomes including mito
        '''
        for rec in self.fasta:
            print(rec.description)

    def print_GC_information(self, window_size = 5000, is_print = True):
        chrom_GC = defaultdict(list)
        chrom_size = defaultdict(int)
        for rec in self.fasta:
            if rec.description.find('SO=chromosome') == -1:
                continue
            chrom = rec.description.split(' | ')[0].split('_')[1]
            length = int(rec.description.split(' | ')[-2].split('=')[1]) 
            for start in range(0, length, window_size):
                end = start + window_size if start + window_size < length else length
                seq = rec.seq[start: end]
                GC = SeqUtils.GC(seq)
                chrom_GC[chrom].append(GC)
            chrom_size[chrom] = length
        if is_print:
            print('maxGC: ', max(max(x) for x in chrom_GC.values()))
            print('minGC: ', min(min(x) for x in chrom_GC.values()))
        return chrom_GC, chrom_size

if __name__ == '__main__':
    genome = General_Genome()
    # genome.print_all_description()
    genome.print_GC_information()

