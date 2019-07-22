library(DESeq2)
library(limma)
library(plyr)

targe.cluster = 1
member = read.csv('D:/Codes/DataLists/Output/member.csv')

d = join_all(lapply(c(1, 4, 5, 7), function(x) {
    tmp = read.table(paste('D:/Codes/DataLists/Output/SRR604939', x, '.dge.txt', sep = ''), header = T)
    colnames(tmp)[-1] = paste(x, 1:(ncol(tmp) - 1), sep = '_')
    tmp
}),
by = 'GENE', type = 'inner') # 15359 * 11269
ieg = as.vector(read.table('D:/Codes/DataLists/Output/IEG.txt')$V1)
rRNA = as.vector(unique(read.table('D:/Codes/DataLists/Output/rRNA.txt')$V1))
d = d[!(d$GENE %in% ieg) & !grepl('mt-', d$GENE) & !(d$GENE %in% rRNA),]
rownames(d) = d$GENE
d$GENE = NULL
d = d[,colSums(d) > 700 & colSums(d != 0) > 450]
d = log(sweep(d, 2, colSums(d), '/') * median(colSums(d)) + 1)

group.list = ifelse(member == targe.cluster, 'target', 'other')
colData = data.frame(row.names = colnames(d), group.list = group.list)
dds = DESeqDataSetFromMatrix(countData = d,
                             colData = colData,
                             design = ~ group.list)
dds = DESeq(dds)
res = results(dds, 
              contrast = c('group.list', 'target', 'other'))
resOrdered = res[order(res$padj),]
nrDEG = as.data.frame(resOrdered)
nrDEG = na.omit(nrDEG)

plot(nrDEG$log2FoldChange, -log10(nrDEG$pvalue))
