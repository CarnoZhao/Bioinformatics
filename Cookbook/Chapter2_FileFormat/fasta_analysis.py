from Bio import Entrez, SeqIO, Seq
from Bio.Alphabet import IUPAC

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
    
    def print_basic_information(self, fasta = '../../../DataLists/CookbookData/NM_002299.fasta'):
        recs = SeqIO.parse(fasta, 'fasta')
        for rec in recs:
            seq = rec.seq
            print(rec.description)
            print(seq[:10])
            print(seq.alphabet)
        return seq

    def transcribe_translate(self, seq):
        '''
        I don't know why need we to change the alphabet from single letter alphabet to unambiguous_dna alphabet
        The original seq can be transcribed/translated too.
        '''
        seq = Seq.Seq(str(seq), IUPAC.unambiguous_dna)
        rna = seq.transcribe()
        protein = seq.translate()
        print("RNA: ", rna)
        print("Protein: ", protein)

if __name__ == '__main__':
    fa = Fasta_Analysis()
    # fia.write_fasta()
    seq = fa.print_basic_information()
    fa.transcribe_translate(seq)
    
