library(dplyr)
library(Seurat)

path = '/mnt/d/Codes/DataLists/GeneMatries/LiverMatrix/'
# path = '/p200/liujiang_group/yinyao/Dataset/Seurat/GeneMatries/LiverMatrix/'

### Process of Raw Data
data = Read10X(data.dir = path)
liver = CreateSeuratObject(counts = data, project = 'Liver', min.cells = 3, min.features = 200)
liver = subset(liver, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
liver = NormalizeData(liver, normalization.method = 'LogNormalize', scale.factor = 10000)
liver = FindVariableFeatures(liver, selection.method = 'vst', nfeatures = 2000)
liver = ScaleData(liver, features = rownames(liver))

### PCA
liver = RunPCA(liver, features = VariableFeatures(object = liver))
# print(liver[["pca"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(liver, dims = 1:2, reduction = "pca")
DimPlot(liver, reduction = "pca")
# DimHeatmap(liver, dims = 1, cells = 500, balanced = T)
# DimHeatmap(liver, dims = 1:15, cells = 500, balanced = TRUE)

### Dimention Determination
liver = JackStraw(liver, num.replicate = 100)
liver = ScoreJackStraw(liver, dims = 1:20)
JackStrawPlot(liver, dims = 1:20)
ElbowPlot(liver)

### Clustering
liver = FindNeighbors(liver, dims = 1:10)
liver = FindClusters(liver, resolution = 0.5)
# head(Idents(liver), 5)

### UMAP
liver = RunUMAP(liver, dims = 1:10)
DimPlot(liver, reduction = "umap")

### tSNE
liver = RunTSNE(liver, dims = 1:10)
DimPlot(liver, reduction = 'tsne')

### Find Markers
liver.markers = FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top1 = liver.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

FeaturePlot(liver, features = top1$gene)


top10 = liver.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(liver, features = top10$gene) # + NoLegend()

new.cluster.ids = 0:13
names(new.cluster.ids) = levels(liver)
liver = RenameIdents(liver, new.cluster.ids)
DimPlot(liver, reduction = "tsne", label = TRUE, pt.size = 0.5, label.size = 7) + NoLegend()
