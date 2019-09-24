import vcf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from variant_calling_format import Variant_Calling_Format
import numpy as np
from collections import defaultdict
import seaborn as sns

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
            if not rec.is_snp or len(rec.ALT) > 1:
                continue
            which_bin = int((rec.POS - start) / total_length * bin_num)
            bins[which_bin] += 1
        return bins

    def bins_MQ0_distribution(self, bin_num = 40, total_length = 200000):
        '''
        plot the median and 75 quantile of MQ0 value for each window, wich window size 2000 out of 200000 total length
        '''
        bins = [[] for _ in range(bin_num)]
        start = None
        for rec in self.v:
            if start == None:
                start = rec.POS
            if not rec.is_snp or len(rec.ALT) > 1:
                continue
            which_bin = int((rec.POS - start) / total_length * bin_num)
            for sample in rec.samples:
                if sample['MQ0'] == None:
                    continue
                bins[which_bin].append(int(sample['MQ0']))
        return bins

    def bins_het_DP_distribution(self):
        '''
        plot the relationsip between hetrozygosity% and DP (read depth), and the relationsip between number of calls and DP
        '''
        het_dic = defaultdict(int)
        total_dic = defaultdict(int)
        for rec in self.v:
            if not rec.is_snp:
                continue
            for sample in rec.samples:
                try:
                    if sample['DP'] == None:
                        continue
                    else:
                        het_dic[int(sample['DP'])] += sample.is_het
                        total_dic[int(sample['DP'])] += 1
                except:
                    pass
        return het_dic, total_dic

    def boxes_EFF_DP(self):
        effs = ['INTERGENIC', 'INTRON', 'NON_SYNONYMOUS_CODING', 'SYNONYMOUS_CODING', 'OTHER']
        dic = dict((eff, []) for eff in effs)
        for rec in self.v:
            if not rec.is_snp:
                continue
            try:
                eff = rec.INFO['EFF'][0].split('(')[0]
                dic[eff].append(int(rec.INFO['DP']))
            except:
                dic['OTHER'].append(int(rec.INFO['DP']))
        return dic

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

    def plot_bins_2(self, name_list = ['centro', 'standard'], bin_num = 40, total_length = 200000):
        bin_list = {}
        for name in name_list:
            vda = Vcf_Data_Analysis('../../../DataLists/CookbookData/%s.vcf.gz' % name)
            bins = vda.bins_MQ0_distribution(bin_num, total_length)
            bin_list[name + ' median'] = [np.median(x) if x else None for x in bins]
            bin_list[name + ' 75 percentile'] = [np.percentile(x, 75) if x else None for x in bins]
        fig, ax = plt.subplots(figsize = (16, 9))
        x_axis = [x * total_length / bin_num for x in range(len(bins))]
        for name, bins in bin_list.items():
            ax.plot(x_axis, bins, 
                    '--' if '75' in name else '-', 
                    label = name, 
                    color = 'b' if name_list[0] in name else 'g'
                    )
        ax.legend()
        ax.set_xlabel('Genomic Location in the Downloaded Segment')
        ax.set_ylabel('MQ0')
        fig.suptitle('Distribution of MQ0 along the Genome', fontsize = 'xx-large')
        plt.savefig('MQ0_distribution.pdf')

    def plot_bins_3(self, name_list = ['centro', 'standard']):
        fig, ax = plt.subplots(figsize = (16, 9))
        ax2 = ax.twinx()
        colors = ['g', 'b']
        for i, name in enumerate(name_list):
            vda = Vcf_Data_Analysis('../../../DataLists/CookbookData/%s.vcf.gz' % name)
            het, total = vda.bins_het_DP_distribution()
            dps = sorted(total.keys())
            ax.plot(dps, [het[dp] / total[dp] for dp in dps], color = colors[i], label = name)
            ax2.plot(dps, [total[dp] for dp in dps], '--', color = colors[i])
        ax.set_xlabel('Sample Read Depth (DP)')
        ax.set_ylabel('Fraction of Heterozygate Calls')
        ax2.set_ylabel('Number of Calls')
        ax.set_xlim(0, 75)
        ax.set_ylim(0, 0.2)
        ax.legend()
        fig.suptitle('Number of Calls per Depth and Fraction of Calls which are Hz', fontsize = 'xx-large')
        plt.savefig('Het_Depth_dsitribution.pdf')
            
    def plot_boxes(self, name = '../../../DataLists/CookbookData/standard.vcf.gz'):
        vda = Vcf_Data_Analysis(name)
        boxes = vda.boxes_EFF_DP()
        fig, ax = plt.subplots(figsize = (16, 9))
        sns.boxplot(data = list(boxes.values()), sym = '', ax = ax)
        ax.set_xticklabels(boxes.keys())
        ax.set_ylabel('DP (Variant)')
        fig.suptitle('Distribution of Variant DP per SNP Site', fontsize = 'xx-large')
        plt.savefig('EFF_boxplot.pdf')

if __name__ == '__main__':
    f = Functions()
    f.plot_boxes()

