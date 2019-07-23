#!/p200/liujiang_group/yinyao/Anaconda/bin/Rscript
#PBS -q bioque
#PBS -l mem=15gb,nodes=1:ppn=5,walltime=20:00:00
#PBS -e /p200/liujiang_group/yinyao/Dataset/Seurat/Logs/Tree.e
#PBS -o /p200/liujiang_group/yinyao/Dataset/Seurat/Logs/Tree.o
#HSCHED -s 1+1+human

library(dplyr)
library(Seurat)

pdf('/p200/liujiang_group/yinyao/Dataset/Seurat/Output.pdf')

### Load Data
path = '/p200/liujiang_group/yinyao/Dataset/Seurat/GeneMatries/LiverMatrix/'
data = Read10X(data.dir = path)
orig.liver = CreateSeuratObject(counts = data, project = 'Orig.Liver', min.cells = 3, min.features = 200)
orig.liver = subset(orig.liver, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

matrix_processing = function(input) {
    liver = NormalizeData(input, normalization.method = 'LogNormalize', scale.factor = 10000)
    liver = FindVariableFeatures(liver, selection.method = 'vst', nfeatures = 2000)
    liver = ScaleData(liver, features = rownames(liver))
    liver = RunPCA(liver, features = VariableFeatures(object = liver))
    liver = FindNeighbors(liver, dims = 1:10)
    liver = FindClusters(liver, resolution = 0.5)
    return(liver)
}


### First Clustering
one_step_cluster = function(input) {
    liver = matrix_processing(input)
    liver.markers = FindAllMarkers(liver, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
    top10 = liver.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(liver, features = top10$gene) + NoLegend()
    fc = top10 %>% group_by(cluster) %>% summarise(mean_avg_logFC = mean(avg_logFC))
    out.cluster = fc$cluster[fc$mean_avg_logFC == max(fc$mean_avg_logFC)]
    return(names(liver@active.ident)[liver@active.ident == out.cluster])
}

liver = orig.liver
clst.array = c()
clst.num = 0
for (i in 0:13) {
    filter = one_step_cluster(liver)
    filter.array = rep(clst.num, length(filter))
    names(filter.array) = filter
    clst.array = c(clst.array, filter.array)
    clst.num = clst.num + 1
    liver = subset(liver, cells = filter, invert = T)
}

orig.liver = matrix_processing(orig.liver)
clusters = clst.array[names(orig.liver@active.ident)]
clusters[is.na(names(clusters))] = clst.num
cells = names(orig.liver@active.ident)
orig.liver@active.ident = factor(clusters)
names(orig.liver@active.ident) = cells
DimPlot(orig.liver)

dev.off()


