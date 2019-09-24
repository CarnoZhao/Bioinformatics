from Bio import SeqIO
from general_genome import General_Genome
import gzip
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

class Low_Quality_Genome():
    '''
    deal with genome fasta file with low quality
    '''
    
    def __init__(self, filename = '../../../DataLists/CookbookData/gambiae.fa.gz'):
        self.gg = General_Genome(gzip.open(filename, 'rt', encoding = 'utf-8'))

    def print_N_information(self):
        '''
        count the number of N fragments in a chromosome and the max length of one N fragment
        '''
        chrom_N = defaultdict(list)
        chrom_size = defaultdict(int)
        for rec in self.gg.fasta:
            chrom = rec.description.split(':')[2]
            if chrom in ('UNKN', 'Y_unplaced'):
                continue
            cnt = 0
            for x in rec.seq:
                if x in ('N', 'n'):
                    cnt += 1
                else:
                    chrom_N[chrom].append(cnt)
                    cnt = 0
            if cnt != 0:
                chrom_N[chrom].append(cnt)
            chrom_size[chrom] = len(rec.seq)
        for chrom in chrom_N:
            print('%s (%s): %%Ns (%.1f), num Ns: %d, max N: %d' % (chrom, chrom_size[chrom], 100 * sum(chrom_N[chrom]) / chrom_size[chrom], len(chrom_N[chrom]), max(chrom_N[chrom])))

    def plot_N_distribution(self):
        '''
        plot the relationship between contig size and fragments of Ns
        '''
        contig_size = []
        contig_N = []
        for rec in self.gg.fasta:
            contig_size.append(len(rec.seq))
            cnt = 0
            for x in rec.seq:
                if x in 'Nn':
                    cnt += 1
            contig_N.append(cnt * 100 / len(rec.seq))
        print('Num of ontigs: ', len(contig_size))
        print('Min size: ', min(contig_size))
        print('Percentile 10: ', np.percentile(contig_size, 10))
        print('Median size: ', np.median(contig_size))
        print('Mean size: ', np.mean(contig_size))
        print('Percentile 90: ', np.percentile(contig_size, 90))
        print('Max size: ', max(contig_size))
        bottom_size = 4800
        top_size = 540000
        fig, ax = plt.subplots(1, 3, figsize = (16, 9), sharey = True, squeeze = False)
        xs, ys = zip(*[(x, y) for x, y in zip(contig_size, contig_N) if x <= bottom_size])
        ax[0, 0].plot(xs, ys, '.')
        xs, ys = zip(*[(x, y) for x, y in zip(contig_size, contig_N) if bottom_size < x <= top_size])
        ax[0, 1].plot(xs, ys, '.')
        ax[0, 1].set_xlim(bottom_size, top_size)
        xs, ys = zip(*[(x, y) for x, y in zip(contig_size, contig_N) if x > top_size])
        ax[0, 2].plot(xs, ys, '.')
        ax[0, 0].set_ylabel('Fraction of Ns', fontsize = 12)
        ax[0, 1].set_xlabel('Contig Size', fontsize = 12)
        fig.suptitle('Fraction of Ns per contig size', fontsize = 26)
        plt.savefig('./Plots/N_contig_size.pdf')


if __name__ == '__main__':
    # lqg = Low_Quality_Genome()
    # lqg.print_N_distribution()
    lqg = Low_Quality_Genome('../../../DataLists/CookbookData/atroparvus.fa.gz')
    lqg.plot_N_distribution()

