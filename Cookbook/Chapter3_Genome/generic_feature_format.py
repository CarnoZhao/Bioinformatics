import gffutils
import sqlite3
from collections import defaultdict
from Bio import Alphabet, Seq, SeqIO
import gzip
from general_genome import General_Genome

class Generic_Feature_Format():
    '''
    use gffutils package to deal with the gff files
    '''

    def __init__(self, filename = '../../../DataLists/CookbookData/gambiae.gff3.gz'):
        try:
            self.gff = gffutils.create_db(filename, '/'.join(filename.split('/')[:-1]) + '/' + filename.split('/')[-1].split('.')[0] + '.db')
        except:
            self.gff = gffutils.FeatureDB('/'.join(filename.split('/')[:-1]) + '/' + filename.split('/')[-1].split('.')[0] + '.db')

    def print_feature_types(self):
        '''
        print the feature type of gff file, such as gene, exon, intron, etc.
        '''
        for feature in self.gff.featuretypes():
            print(feature, self.gff.count_features_of_type(feature))

    def print_all_contigs(self):
        '''
        print all the contigs, that is some well-built chromosomes (including mito) and some unknown parts
        '''
        for contig in self.gff.features_of_type('contig'):
            print(contig)

    def deeper_analysis(self):
        '''
        1. print the number of genes per contig
        2. print the gene with most exons
        3. print the max length (or other length statistics) of genes
        4. print the distribution of the number of mRNAs (all)
        5. print the distribution of the number of exons (all)
        '''
        num_mRNAs = defaultdict(int)
        num_exons = defaultdict(int)
        max_span = 0
        max_exon = 0
        for contig in self.gff.features_of_type('contig'):
            cnt = 0
            for gene in self.gff.region((contig.id, contig.start, contig.end), featuretype = 'gene'):
                cnt += 1
                span = abs(gene.start - gene.end)
                if span > max_span:
                    max_span = span
                    max_span_gene = gene
                mRNAs = list(self.gff.children(gene, featuretype = 'mRNA'))
                num_mRNAs[len(mRNAs)] += 1
                checks = mRNAs if mRNAs else [gene]
                for check in checks:
                    exons = list(self.gff.children(check, featuretype = 'exon'))
                    num_exons[len(exons)] += 1
                    if len(exons) > max_exon:
                        max_exon = len(exons)
                        max_exon_gene = gene
            print('contig %s: %d genes' % (contig.seqid, cnt))
        print('Max number of exons: %s (%d)' % (max_exon_gene.id, max_exon))
        print('Max span: %s (%d)' % (max_span_gene.id, max_span))
        print('Number of mRNAs: %s' % num_mRNAs)
        print('Number of exons: %s' % num_exons)

    def extract_genes(self, gene_id = 'AGAP004707'):
        gene = self.gff[gene_id]
        print('gene: ', gene)
        print('gene seqid and strand: ', gene.seqid, gene.strand)
        genome = General_Genome(gzip.open('../../../DataLists/CookbookData/gambiae.fa.gz', 'rt', encoding = 'utf-8'))
        for rec in genome.fasta:
            print('record description: ', rec.description)
            if rec.description.split(':')[2] == gene.seqid:
                seq = rec.seq
                break
        print(seq.alphabet)
        mRNAs = self.gff.children(gene, featuretype = 'mRNA')
        for mRNA in mRNAs:
            print(mRNA.id)
            if mRNA.id.endswith('RA'):
                break
        CDSs = self.gff.children(mRNA, featuretype = 'CDS', order_by = 'start')
        gene_seq = self.get_sequence(seq, CDSs, gene.strand)
        print(len(gene_seq), gene_seq[:10] + '...' + gene_seq[-10:])
        prot = gene_seq.translate()
        print(len(prot), prot[:10] + '...' + prot[-10:])

    def get_sequence(self, seq, CDSs, strand):
        retseq = Seq.Seq('', alphabet = Alphabet.IUPAC.unambiguous_dna)
        for cds in CDSs:
            inner = Seq.Seq(str(seq[cds.start - 1: cds.end]), alphabet = Alphabet.IUPAC.unambiguous_dna)
            retseq += inner
        return retseq if strand == '+' else retseq.reverse_complement()

if __name__ == '__main__':
    gff = Generic_Feature_Format()
    # gff.print_feature_types()
    # gff.print_all_contigs()
    # gff.deeper_analysis()
    gff.extract_genes()
