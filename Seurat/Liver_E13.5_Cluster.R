library(dplyr)
library(Seurat)
library(ggplot2)

path = '/mnt/d/Codes/DataLists/GeneMatries/LiverMatrix/'
data = Read10X(data.dir = path)
liver = CreateSeuratObject(counts = data, project = 'Liver', min.cells = 3, min.features = 200)
liver = subset(liver, subset = nFeature_RNA > 200)

matrix_processing = function(input, reso) {
  liver = NormalizeData(input, normalization.method = 'LogNormalize', scale.factor = 10000)
  liver = FindVariableFeatures(liver, selection.method = 'vst', nfeatures = 2000)
  liver = ScaleData(liver, features = rownames(liver))
  liver = RunPCA(liver, features = VariableFeatures(object = liver))
  liver = RunTSNE(liver, dims = 1:10)
  liver = RunUMAP(liver, dims = 1:10)
  liver = FindNeighbors(liver, dims = 1:10)
  liver = FindClusters(liver, resolution = reso)
  return(liver)
}

### first clustering
liver = matrix_processing(liver, 0.4)


liver.markers = FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top12 = liver.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_logFC)
graph = DoHeatmap(liver, features = top12$gene) + NoLegend()
plot(graph)

DimPlot(liver, reduction = 'umap', label = T)
