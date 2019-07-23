#!/p200/liujiang_group/yinyao/Anaconda/bin/Rscript
#PBS -q bioque
#PBS -l mem=15gb,nodes=1:ppn=5,walltime=20:00:00
#PBS -e /p200/liujiang_group/yinyao/Dataset/Seurat/Logs/Tree.e
#PBS -o /p200/liujiang_group/yinyao/Dataset/Seurat/Logs/Tree.o
#HSCHED -s 1+1+human

library(dplyr)
library(Seurat)

# pdf('/p200/liujiang_group/yinyao/Dataset/Seurat/Output.pdf')
pdf('/mnt/d/Codes/DataLists/GeneMatries/LiverOutput/Output.pdf')

### Load Data
# path = '/p200/liujiang_group/yinyao/Dataset/Seurat/GeneMatries/LiverMatrix/'
path = '/mnt/d/Codes/DataLists/GeneMatries/LiverMatrix/'
data = Read10X(data.dir = path)
orig.liver = CreateSeuratObject(counts = data, project = 'Orig.Liver', min.cells = 3, min.features = 200)
orig.liver = subset(orig.liver, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

matrix_processing = function(input) {
    liver = NormalizeData(input, normalization.method = 'LogNormalize', scale.factor = 10000)
    liver = FindVariableFeatures(liver, selection.method = 'vst', nfeatures = 2000)
    liver = ScaleData(liver, features = rownames(liver))
    liver = RunPCA(liver, features = VariableFeatures(object = liver))
    liver = RunTSNE(liver, dims = 1:10)
    liver = RunUMAP(liver, dims = 1:10)
    liver = FindNeighbors(liver, dims = 1:10)
    liver = FindClusters(liver, resolution = 0.5)
    return(liver)
}


### First Clustering
one_step_cluster = function(input, i) {
    liver = matrix_processing(input)
    liver.markers = FindAllMarkers(liver, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
    top10 = liver.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(liver, features = top10$gene) + NoLegend()
    write.csv(top10, file = paste('/mnt/d/Codes/DataLists/GeneMatries/LiverOutput/Liver', i, '.csv', sep = ''), row.names = T)
    fc = top10 %>% group_by(cluster) %>% summarise(mean_avg_logFC = mean(avg_logFC))
    out.cluster = fc$cluster[rev(order(fc$mean_avg_logFC))[1:2]]
    saveRDS(liver, paste('/mnt/d/Codes/DataLists/GeneMatries/LiverOutput/Liver', i, '.rds', sep = ''))
    return(liver@active.ident[liver@active.ident %in% out.cluster])
}

liver = orig.liver
ident = orig.liver@active.ident
i = 0
filter = NULL
while (length(liver@active.ident) != 0) {
    filter = one_step_cluster(liver, i)
    ident[names(filter)] = filter + i
    liver = subset(liver, cells = filter, invert = T)
    i = i + 2
}

orig.liver = matrix_processing(orig.liver)
orig.liver@active.ident = as.factor(ident)
DimPlot(orig.liver, reduction = 'tsne')

dev.off()


