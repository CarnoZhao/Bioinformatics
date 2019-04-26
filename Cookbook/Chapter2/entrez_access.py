from Bio import Entrez, SeqIO

class Entrez_Access():
    '''
    Accessing the Entrez and other NCBI database using Biopython
    '''

    def __init__(self):
        Entrez.email = 'xz2827@columbia.edu'

    def search_gene(self, gene = 'TP53', organism = 'Homo Sapiens'):
        '''
        Give gene name and organism name
        Return the results of entrez search
        '''
        handle = Entrez.esearch(db = 'nucleotide', term = '%s[GENE] AND "%s"[ORGANISM]' % (gene, organism))
        rec_list = Entrez.read(handle)
        if rec_list['RetMax'] < rec_list['Count']:
            handle = Entrez.esearch(db = 'nucleotide', term = '%s[GENE] AND "%s"[ORGANISM]' % (gene, organism), retmax = rec_list['Count'])
            rec_list = Entrez.read(handle)
        return rec_list

    def parse_esearch_result(self, rec_list):
        '''
        give the rec_list of `search_gene`
        return certain result
        '''
        id_list = rec_list['IdList']
        hdl = Entrez.efetch(db = 'nucleotide', id = id_list, rettype = 'gb') # genbank
        recs = list(SeqIO.parse(hdl, 'gb'))
        rec = recs[0]
        return rec

    def print_gene_information(self, rec):
        print("name", rec.name)
        print("description", rec.description)
        print('Features:')
        for feature in rec.features:
            if feature.type == 'gene':
                print('gene:', feature.qualifiers['gene'])
            elif feature.type == 'exon':
                loc = feature.location
                print('exon:', loc.start, loc.end, loc.strand)
            else:
                print('skip %s' % feature.type)

    def print_gene_annotation(self, rec):
        for name, value in rec.annotations.items():
            print('%s == %s' % (name, value))

    def print_seq(self, rec):
        print(len(rec.seq))
        print(rec.seq)


if __name__ == '__main__':
    ea = Entrez_Access()
    rec_list = ea.search_gene()
    rec = ea.parse_esearch_result(rec_list)
    # ea.print_gene_information(rec)
    # ea.print_gene_annotation(rec)
    ea.print_seq(rec)

