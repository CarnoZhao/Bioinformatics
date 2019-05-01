from ensembl_rest_api import Ensembl_Rest_Api

class Gene_Ontology(Ensembl_Rest_Api):
	'''
	use pygraphviz and ensembl api to visualize the gene ontology 
	'''

	def __init__(self, gene_id = "ENSG00000115850"): # lactase gene
		super().__init__()
		self.gene_id = gene_id
		self.refs = self.do_request('xrefs/id', gene_id, external_db = 'GO', all_levels = '1')

	def general_information(self):
		print(len(self.refs))
		print(self.refs[0].keys())
		for ref in self.refs:
			go_id = ref['primary_id']
			details = self.do_request('ontology/id', go_id)
			print('%s %s %s' % (go_id, details['namespace'], ref['description']))
			print('%s\n' % details['definition'])

	def detail_about_go_id(self, go_id = 'GO:0000016'):
		'''
		get the detail information of one go id node, including its parent
		'''
		data = self.do_request('ontology/id', go_id)
		for k, v in data.items():
			if k == 'parents':
				for parent in v:
					print(parent)
					parent_id = parent['accession']
			else:
				print('%s: %s' % (k ,str(v)))
		parent_data = self.do_request('ontology/id', parent_id)
		print(parent_id, len(parent_data['children']))

	def get_all_ancestors(self, go_id = 'GO:0000016'):
		'''
		get all the ancestors of a certain go id, stop until 3 roots (Molecular funtion, cell component, biological processes)
		'''
		refs = self.do_request('ontology/ancestors/chart', go_id)
		for go, entry in refs.items():
			print(go)
			term = entry['term']
			print('%s: %s' % (term['name'], term['definition']))
			is_a = entry.get('is_a', [])
			print('\t is a: %s\n' % ', '.join([x['accession'] for x in is_a]))

	def get_upper(self, go_id):
		

if __name__ == '__main__':
	go = Gene_Ontology()
	# go.general_information()
	# go.detail_about_go_id()
	go.get_all_ancestors()