library(gmodels)
library(ggplot2)
library(plyr)

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
# pca = fast.prcomp(t(d), center = FALSE)
# rm(pca)
PC1_50 = as.data.frame(pca$x[,1:50])
write.csv(PC1_50, file = 'D:/Codes/DataLists/Output/PC1_50.csv', row.names = T)
# PC1_50$member = read.csv('D:/Codes/DataLists/Output/member.csv')
ggplot(PC1_50, aes(x = PC1, y = PC2)) + geom_point()
