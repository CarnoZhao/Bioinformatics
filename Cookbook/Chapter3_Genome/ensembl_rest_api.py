import requests

class Ensembl_Rest_Api():
    '''
    using web requests package to get data from ensembl database, here only search for orthologues
    '''
    
    def __init__(self):
        self.url = 'http://rest.ensembl.org'

    def do_request(self, service, *args, **kwargs):
        url_params = ''
        for a in args:
            if a is not None:
                url_params += '/' + a
        req = requests.get('%s/%s%s' % (self.url, service, url_params), params = kwargs, headers = {'Content-Type': 'application/json'})
        if not req.ok:
            req.raise_for_status()
        return req.json()

    def get_species(self):
        species = self.do_request('info/species')
        for i, specie in enumerate(species['species']):
            print(i, specie['name'])

    def get_human_database(self):
        answer = self.do_request('info/external_dbs', 'homo_sapiens', filter = 'HGNC%')
        print(answer)

    def get_gene_id(self, organism = 'homo_sapiens', gene_name = 'LCT', is_print = True):
        answer = self.do_request('lookup/symbol', organism, gene_name)
        if is_print:
            print(answer)
        return answer['id']

    def get_gene_seq(self, organism = 'homo_sapiens', gene_name = 'LCT'):
        answer = self.do_request('sequence/id', self.get_gene_id(is_print = False))
        print(answer)
        # extra reference database
        xrefs = self.do_request('xrefs/id', self.get_gene_id(is_print = False))
        for xref in xrefs:
            print(xref['db_display_name'])
            print(xref)

    def get_orthologues(self, org1 = 'homo_sapiens', org2 = 'equus_caballus', gene_name = 'LCT'):
        hom_response = self.do_request('homology/id', self.get_gene_id(organism = org1, gene_name = gene_name, is_print = False), type = 'orthologues', sequence = 'none')
        homologies = hom_response['data'][0]['homologies']
        for homology in homologies:
            print(homology['target']['species'])
            if homology['target']['species'] != org2:
                continue
            print(homology)
            print(homology['taxonomy_level'])
            horse_id = homology['target']['id']
        horse_req = self.do_request('lookup/id', horse_id)
        print(horse_req)

if __name__ == '__main__':
    era = Ensembl_Rest_Api()
    # era.get_species()
    # era.get_human_database()
    # era.get_gene_id()
    # era.get_gene_seq()
    era.get_orthologues()
