library(dplyr)
library(Seurat)
library(ggplot2)
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
load.orig = function(){
    path = '/mnt/d/Codes/DataLists/GeneMatries/LiverMatrix/'
    data = Read10X(data.dir = path)
    liver = CreateSeuratObject(counts = data, project = 'Liver', min.cells = 3, min.features = 200)
    liver[["percent.mt"]] = PercentageFeatureSet(liver, pattern = "^MT-")
    liver = subset(liver, subset = nFeature_RNA > 500 & percent.mt < 10)
    liver
}

liver = matrix_processing(liver, 0.8)

saveRDS(liver, '/mnt/d/Codes/DataLists/GeneMatries/Liver.rds', compress = T)

liver = readRDS('/mnt/d/Codes/DataLists/GeneMatries/Liver.rds')
markers = read.csv('/mnt/d/Codes/DataLists/GeneMatries/markers.csv', header = T, row.names = 1)
markers$cluster = as.factor(markers$cluster)
pvalue.markers = markers[markers$p_val_adj < 0.05,]
um = DimPlot(liver, reduction = 'umap', label = T, label.size = 5)
ts = DimPlot(liver, reduction = 'tsne', label = T, label.size = 5)

markers = FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10$gene = as.character(top10$gene)
heat = DoHeatmap(liver, features = top10$gene) + NoLegend()

map = c('6'  = 'Hepatoblast',
        '13' = 'Hepatocyte',
        '12' = 'Macrophage',
        '14' = 'Megakaryocyte',
        '0'  = 'Erythroid',
        '1'  = 'Erythroid',
        '5'  = 'Erythroid',
        '2'  = 'Erythroid Progenitor',
        '4'  = 'Erythroid Progenitor',
        '8'  = 'Erythroid Progenitor',
        '9'  = 'Erythroid Progenitor',
        '7'  = 'Erythroid Progenitor')

g = function(i){
    gs = pvalue.markers[pvalue.markers$cluster %in% i,]
    ge = pvalue.markers[!pvalue.markers$cluster %in% i,]
    View(gs[!gs$gene %in% ge$gene,])
}

### 11
liver = load.orig()
liver = subset(liver, cells = cell11)
liver = matrix_processing(liver, 0.6)
markers = FindAllMarkers(liver, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
cell_Endo = names(liver@active.ident[liver$seurat_clusters == 0])
cell_Mesen = names(liver@active.ident[liver$seurat_clusters == 1])

### 10
liver = load.orig()
liver = subset(liver, cells = cell10)
liver = matrix_processing(liver, 0.6)
markers = FindAllMarkers(liver, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
top10 = markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(liver, features = top10$gene)
cell_Neutr = names(liver@active.ident[liver$seurat_clusters == 2])
cell_Stem = names(liver@active.ident[liver$seurat_clusters == 0])
cell_EryPro = names(liver@active.ident[liver$seurat_clusters == 1])


back.up = liver@active.ident
ident = as.character(back.up)
ident = ifelse(ident %in% names(map), map[ident], ident)
ident = ifelse(names(back.up) %in% cell_Endo, 'Endothelial', ident)
ident = ifelse(names(back.up) %in% cell_EryPro, 'Erythroid Progenitor', ident)
ident = ifelse(names(back.up) %in% cell_Mesen, 'Mesenchymal', ident)
ident = ifelse(names(back.up) %in% cell_Neutr, 'Neutrophil', ident)
ident = ifelse(names(back.up) %in% cell_Stem, 'Hematopoietic', ident)
ident = as.factor(ident)
names(ident) = names(back.up)
liver@active.ident = ident
liver$seurat_clusters = ident
DimPlot(liver, reduction = 'umap', label = T)
DimPlot(subset(liver, subset = seurat_clusters != 3), reduction = 'umap', label = T)


