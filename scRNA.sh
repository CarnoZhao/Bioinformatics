########
# name #
########

species=''

############
# software #
############

picard='' # path to picard.jar



#########
# files #
#########

fasta='' # reference sequence of organism
dict='' # created from the fasta
# created by Picard CreateSequenceDictionary
gtf='' # locatiion of genomic features
refFlat='' # same information but defferent format of gtf file
# created by ConvertToRefFlat in DropSeq
genes_intervals='' # optional
exons_intervals='' # optional
rRNA_intervals='' # optional
reduced_gtf='' # optional

# mouse data available in ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz

# human data available in http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1629193

###############
# preparation #
###############

# create seq-dict
java -jar $picard CreateSequenceDictionary \
    REFERENCE=$fasta \
    OUTPUT=$dict \
    SPECIES=$species

# create RefFlat
