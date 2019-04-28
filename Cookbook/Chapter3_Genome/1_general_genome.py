from Bio import SeqIO, SeqUtils
from collections import defaultdict
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome

class General_Genome():
    '''
    general information about a genome fasta file
    '''
    def __init__(self, filename = '../../../DataLists/CookbookData/PlasmoDB-9.3_Pfalciparum3D7_Genome.fasta'):
        self.fasta = SeqIO.parse(filename, 'fasta')
        
    def print_all_description(self):
        '''
        print the description of all fasta records, which is separated to different chromosomes including mito
        '''
        for rec in self.fasta:
            print(rec.description)

    def print_GC_information(self, window_size = 50000, is_print = True):
        chrom_GC = defaultdict(list)
        chrom_size = defaultdict(int)
        for rec in self.fasta:
            if rec.description.find('SO=chromosome') == -1:
                continue
            chrom = int(rec.description.split(' | ')[0].split('_')[1])
            length = int(rec.description.split(' | ')[-2].split('=')[1]) 
            for start in range(0, length, window_size):
                end = start + window_size if start + window_size < length else length
                seq = rec.seq[start: end]
                GC = SeqUtils.GC(seq)
                chrom_GC[chrom].append(GC)
            chrom_size[chrom] = length
        if is_print:
            print('maxGC: ', max(max(x) for x in chrom_GC.values()))
            print('minGC: ', min(min(x) for x in chrom_GC.values()))
        return chrom_GC, chrom_size

    def plot_chromosomes_GC(self, window_size = 50000):
        chrom_GC, chrom_size = self.print_GC_information(window_size = window_size, is_print = False)
        chroms = list(chrom_GC.keys())
        chroms.sort()
        longest_chrom = max(chrom_size.values())
        fig = BasicChromosome.Organism(output_format = 'pdf')
        fig.page_size = (29.7 * cm, 21 * cm)
        telomere_length = 10
        bottom_GC = 17.5
        top_GC = 22
        for chrom in chroms:
            size = chrom_size[chrom]
            representation = BasicChromosome.Chromosome('Cr %d' % chrom)
            representation.scale_num = longest_chrom
            tel = BasicChromosome.TelomereSegment()
            tel.scale = telomere_length
            representation.add(tel)
            num_blocks = len(chrom_GC[chrom])
            for block, gc in enumerate(chrom_GC[chrom]):
                body = BasicChromosome.ChromosomeSegment()
                if gc > top_GC:
                    body.fill_color = colors.Color(1, 0, 0)
                elif gc < bottom_GC:
                    body.fill_color = colors.Color(1, 1, 0)
                else:
                    col = (gc- bottom_GC) / (top_GC - bottom_GC)
                    body.fill_color = colors.Color(col, col, 1)
                if block < len(chrom_GC) - 1:
                    body.scale = window_size
                else:
                    body.scale = size % window_size
                representation.add(body)
            tel = BasicChromosome.TelomereSegment(inverted = True)
            tel.scale = telomere_length
            representation.add(tel)
            fig.add(representation)
        fig.draw('GC_chromosome_plot.pdf', 'Plasimodium falciparum')

if __name__ == '__main__':
    genome = General_Genome()
    # genome.print_all_description()
    # genome.print_GC_information()
    genome.plot_chromosomes_GC()

