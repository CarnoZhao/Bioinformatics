from Bio import SeqIO
from 1_general_genome import General_Genome
import gzip

class Low_Quality_Genome():
    '''
    deal with genome fasta file with low quality
    '''
    
    def __init__(self):
        pass

if __name__ == '__main__':
    gg = General_Genome(gzip.open('../../../DataLists/CookbookData/gambiae.fa.gz', 'rt', encoding = 'utf-8'))
    gg.print_information()
