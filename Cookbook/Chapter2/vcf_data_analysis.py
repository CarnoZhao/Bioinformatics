import vcf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from variant_calling_format import Variant_Calling_Format

class Vcf_Data_Analysis(Variant_Calling_Format):    
    '''
    some data analysis of vcf file, based on information given in `Variant_Calling_Format` class
    '''

    def bins_bi_allelic_distribution(self, bin_num = 100, total_length = 200000):
        '''
        plot the distribution of bi-allelic (is_snp and len(ALT) == 1)
        '''
        bins = [0 for _ in range(bin_num)]
        start = None
        for rec in self.v:
            if start == None:
                start = rec.POS
            if not(rec.is_snp and len(rec.ALT) == 1):
                continue
            which_bin = int((rec.POS - start) / total_length * bin_num)
            bins[which_bin] += 1
        return bins

    def bins_MQ0_distribution(self, bin_num = 100, total_length = 200000):
        '''
        plot the median and 75 quantile of MQ0 value for each window, wich window size 2000 out of 200000 total length
        '''
        bins = [[] for _ in range(bin_num)]
        start = None
        for rec in self.v:
            if start == None:
                start = rec.POS
            if not (rec.is_snp and len(rec.ALT) == 1):
                continue
            which_bin = int((rec.POS - start) / total_length * bin_num)
            for sample in rec.samples:
                if sample['MQ0'] == None:
                    continue
                bins[which_bin].append(int(sample['MQ0'])


class Functions():
    def __init__(self):
        pass

    def plot_bins_1(self, name_list = ['centro', 'standard'], bin_num = 100, total_length = 200000):
        '''
        given a list of bins, with same length, use matplotlib to plot
        '''
        bin_list = []
        for name in name_list:
            vda = Vcf_Data_Analysis('../../../DataLists/CookbookData/%s.vcf.gz' % name)
            bins = vda.bins_bi_allelic_distribution()
            bin_list.append(bins)
        fig, ax = plt.subplots(figsize = (16, 9))
        for bins in bin_list:
            ax.plot([x * total_length / bin_num for x in range(len(bins))], bins)
        ax.legend(labels = [name + '.vcf.gz' for name in name_list])
        ax.set_xlabel('Genomic Location in the Downloaded Segment')
        ax.set_ylabel('Number of Variant Site (bi-allelic SNPs)')
        fig.suptitle('Number of Bi-allelic SNPs along the Genome', fontsize = 'xx-large')
        plt.savefig('bi_allelic_distribution.pdf')

if __name__ == '__main__':
    f = Functions()
    f.plot_bins()

