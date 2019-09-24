library(dplyr)
library(Seurat)

path = '/mnt/d/Codes/DataLists/GeneMatries/PBMC/'
# path = '/p200/liujiang_group/yinyao/Dataset/Seurat/GeneMatries/PBMC/'

data = Read10X(data.dir = path)
pbmc = CreateSeuratObject(counts = data, project = 'pbmc3k', min.cells = 3, min.features = 200)

pbmc[["Percent_mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "Percent_mt"))

plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "Percent_mt")
plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & Percent_mt < 5)

pbmc = NormalizeData(pbmc, normalization.method = 'LogNormalize', scale.factor = 10000)

pbmc = FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)
top10 = head(VariableFeatures(pbmc), 10)
plot1 = VariableFeaturePlot(pbmc)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

pbmc = ScaleData(pbmc)

pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = T)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

pbmc = JackStraw(pbmc, num.replicate = 100)
pbmc = ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(pbmc)

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)

pbmc = RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))


top10 = pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids = c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) = levels(pbmc)
pbmc = RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
