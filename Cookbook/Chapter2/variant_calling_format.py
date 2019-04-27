import vcf
from collections import defaultdict

class Variant_Calling_Format():
    '''
    basic analysis of vcf file using PyVcf (pip install pyvcf, import vcf)
    '''

    def __init__(self, filename = '../../../DataLists/CookbookData/genotypes.vcf.gz'):
        self.v = vcf.Reader(filename = filename)

    def print_information(self):
        print('Variant Level Information')
        infos= self.v.infos
        for info in infos:
            print(info)
        print('Sample Level Information')
        fmts = self.v.formats
        for fmt in fmts:
            print(fmt)

    def print_information_one_record(self):
        rec = next(self.v)
        print(rec.CHROM, rec.POS, rec.ID, rec.REF, rec.ALT, rec.QUAL, rec.FILTER)
        print(rec.INFO)
        print(rec.FORMAT)
        samples = rec.samples
        print(len(samples))
        sample = samples[0]
        print(sample.called, sample.gt_alleles, sample.is_het, sample.phased)
        print(int(sample['DP']))

    def print_variant_type(self):
        my_type = defaultdict(int)
        num_alts = defaultdict(int)
        for rec in self.v:
            my_type[rec.var_type, rec.var_subtype] += 1
            if rec.is_snp:
                num_alts[len(rec.ALT)] += 1
        print(my_type)
        print(num_alts)

    def 

if __name__ == '__main__':
    VCF = Variant_Calling_Format()
    # VCF.print_information()
    # VCF.print_information_one_record()
    VCF.print_variant_type()

