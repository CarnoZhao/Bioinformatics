from Bio import SeqIO
import gzip
from collections import defaultdict, Counter
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns

class Fastq_Analysis():
    '''
    some simple analysis of fastq file, including the qualities and the sequence itself
    '''

    def __init__(self, filename = '../../../DataLists/CookbookData/SRR003265.filt.fastq.gz'):
        self.filename = filename
        self.recs = SeqIO.parse(gzip.open(filename, 'rt', encoding = 'utf-8'), 'fastq')
    
    def print_simple_information(self):
        rec = next(self.recs)
        print('rec id:', rec.id, rec.description, rec.seq)
        print(rec.letter_annotations)

    def count_base_proportion(self):
        '''
        Count the total bases `ACGT` proportion of the whole fastq file
        '''
        cnt = defaultdict(int)
        for rec in self.recs:
            counter = Counter(rec.seq)
            for base in counter:
                cnt[base] += counter[base]
        total = sum(cnt.values())
        for base in cnt:
            print('%s: %.2f%%, %d' % (base, cnt[base] / total, cnt[base]))

    def plot_N_distribution(self):
        '''
        plot the distribution of N base according to the position on the reads of fastq file
        '''
        cnt = defaultdict(int)
        for rec in self.recs:
            for i, base in enumerate(rec.seq):
                if base == 'N':
                    cnt[i] += 1
        poses = list(range(1, max(cnt.keys()) + 2))
        cnt = [cnt[pos - 1] for pos in poses]
        fig, ax = plt.subplots(figsize = (16, 9))
        ax.plot(poses, cnt)
        ax.set_xlim(1, max(poses))
        plt.savefig('N_distribution_fastq.pdf')

    def plot_phred_distribution(self):
        '''
        plot the distribution of phred quality according to the position on the reads of fastq file
        '''
        cnt = defaultdict(list)
        for rec in self.recs:
            for i, p in enumerate(rec.letter_annotations['phred_quality']):
                if i < 25 or p == 40:
                    continue
                cnt[i].append(p)
        cnt = [cnt[i - 1] for i in range(26, max(cnt.keys()) + 2)]
        fid, ax = plt.subplots(figsize = (16, 9))
        sns.boxplot(data = cnt, ax = ax)
        ax.set_xticklabels([str(x) for x in range(26, len(cnt) + 26)])
        plt.savefig('Phred_distribution_fastq.pdf')

if __name__ == '__main__':
    fq = Fastq_Analysis()
    # fq.print_simple_information()
    # fq.count_base_proportion()
    # fq.plot_N_distribution()
    fq.plot_phred_distribution()
    
