########
# name #
########

species=''
filename=''
adapters=''

############
# software #
############

java=''
picard='' # path to picard.jar
dropseq='' # path to dropseq.jar 
star=''

#########
# files #
#########

fasta='' # reference sequence of organism
dict='' # created from the fasta, created by Picard CreateSequenceDictionary
gtf='' # locatiion of genomic features
refFlat='' # same information but defferent format, created by ConvertToRefFlat in DropSeq
genes_intervals='' # optional
exons_intervals='' # optional
rRNA_intervals='' # optional
reduced_gtf='' # optional
star_genome_dir=''
fastq1=''
fastq2=''
unaligned_bam=''
unaligned_tagged_cell_bam=''
unaligned_tagged_cell_mole_bam=''
unaligned_tagged_filtered_bam=''
unaligned_tagged_trimmed_bam=''
unaligned_final_bam=''
unaligned_final_fastq=''
aligned_sam=''
aligned_sorted_bam=''
merged_bam=''
exon_tagged_bam=''
exon_tagged_dge=''
dge_summary=''
tagged_cell_summary=''
tagged_mole_summary=''
adapter_trimmed_summary=''
polyA_trimmed_summary=''



# mouse data available in 
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz
# human data available in 
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1629193

###############
# preparation #
###############

# create unaligned bam
$java -jar $picard FastqToSam \
        FASTQ=$fastq1 \
        FASTQ2=$fastq2 \
        OUTPUT=$unaligned_bam \
        SM=$filename # ???

# create seq-dict
$java -jar $picard CreateSequenceDictionary \
        REFERENCE=$fasta \
        OUTPUT=$dict \
        SPECIES=$species

# create RefFlat
$java -jar $dropseq ConvertToRefFlat \
        ANNOTATIONS_FILE=$gtf \
        SEQUENCE_DICTIONARY=$dict \
        OUPUT=$refFlat

########
# trim #
########

# cell barcode 1-12
$java -jar $dropseq TagBamWithReadSequenceExtended \
        INPUT=$unaligned_bam \
        OUTPUT=$unaligned_tagged_cell_bam \
        SUMMARY=$tagged_cell_summary \
        BASE_RANGE=1-12 \
        BASE_QUALITY=10 \
        BARCODED_READ=1 \
        DISCARD_READ=False \
        TAG_NAME=XC \
        NUM_BASES_BELOW_QUALITY=1

# molecular barcode 13-20
$java -jar $dropseq TagBamWithReadSequenceExtended \
        INPUT=$unaligned_tagged_cell_bam \
        OUTPUT=$unaligned_tagged_cell_mole_bam \
        SUMMARY=$tagged_mole_summary \
        BASE_RANGE=13-20 \
        BASE_QUALITY=10 \
        BARCODED_READ=1 \
        DISCARD_READ=True \
        TAG_NAME=XM \
        NUM_BASES_BELOW_QUALITY=1

# filter bam, to remove low quality tagged when dealing with barcode
$java -jar $dropseq FilterBam \
        TAG_REJECT=XQ \
        INPUT=$unaligned_tagged_cell_mole_bam \
        OUTPUT=$unaligned_tagged_filtered_bam

# trim adapters
$java -jar $dropseq TrimStartingSequence \
        INPUT=$unaligned_tagged_filtered_bam \
        OUTPUT=$unaligned_tagged_trimmed_bam \
        OUTPUT_SUMMARY=$adapter_trimmed_summary \
        SEQUENCE=$adapters \
        MISMATCHES=0 \
        NUM_BASES=5

# trim polyA
$java -jar $dropeq PolyATrimmer \
        INPUT=$unaligned_tagged_trimmed_bam \
        OUTPUT=$unaligned_final_bam \
        OUTPUT_SUMMARY=$polyA_trimmed_summary \
        MISMATCH=0 \
        NUM_BASES=6 \
        USE_NEW_TRIMMER=true

# bam to fastq
$java -jar $picard SamToFastq \
        INPUT=$unaligned_final_bam \
        FASTQ=$unaligned_final_fastq

#############
# alignment #
#############

# aligned to genome
$star --genomeDir $star_genome_dir \
        --readFilesIn $unaligned_final_fastq \
        --outFileNamePrefix $aligned_sam

# sorted by query name
$java -jar $picard SortSam \
        I=$aligned_sam \
        O=$aligned_sorted_bam \
        SO=queryname

# merged mapped and unmapped bam to restore the tags
$java -jar $picard MergeBamAlignment \
        REFERENCE_SEQUENCE=$fasta \
        UNMAPPED_BAM=$unaligned_final_bam \
        ALIGNED_BAM=$aligned_sorted_bam \
        OUTPUT=$merged_bam \
        INCLUDE_SECONDARY_ALIGNMENTS=false \
        PAIRED_RUN=false

###############
# tag elements #
###############


$java -jar $dropseq TagReadWithGeneFunction \
        I=$merged_bam \
        O=$exon_tagged_bam \
        ANNOTATIONS_FILE=$refFlat

#################
# detect errors #
#################

# $java -jar $dropseq DetectBeadSubstitutionErrors \
#         I=$exon_tagged_bam \
#         O=$clean_substitution_bam \
#         OUTPUT_REPORT=$clean_substitution_report

# $java -jar $dropseq DetectBeadSynthesisErrors \
#         I=$clean_substitution_bam \
#         O=$clean_bam \
#         REPORT=$clean_report \
#         OUTPUT_STATS=$clean_stats \
#         SUMMARY=$clean_summary \
#         PRIMER_SEQUENCE=

###############
# extract DGE #
###############

$java -jar $dropseq DigitalExpression \
        I=$exon_tagged_bam \
        O=$exon_tagged_dge \
        SUMMARY=$dge_summary \
        NUM_CORE_BARCODES=100
        # ???