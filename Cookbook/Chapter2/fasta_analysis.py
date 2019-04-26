from Bio import Entrez, SeqIO

class Fasta_Analysis():
    '''
    basic analysis of fasta file using biopython
    '''

    def __init__(self):
        Entrez.email = 'xz2827@columbia.edu'

    def write_fasta(self, accession = ['NM_002299']):
        '''
        given a list of accession number
        write a fasta file in the data dir
        '''
        hdl = Entrez.efetch(db = 'nucleotide', id = accession, rettype = 'fasta')
        seq = SeqIO.read(hdl, 'fasta')
        fasta = open('../../../DataLists/CookbookData/%s.fasta' % accession[0], 'w')
        seq = seq[11:5795]
        SeqIO.write(seq, fasta, 'fasta')
        fasta.close()
    

if __name__ == '__main__':
    fa = Fasta_Analysis()
    fa.write_fasta()
