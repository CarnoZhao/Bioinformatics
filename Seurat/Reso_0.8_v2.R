library(Seurat)
library(ggplot2)
setwd('/mnt/d/Codes/DataLists/GeneMatries/')
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

liver = readRDS('Liver.rds')
d = DimPlot(liver, reduction = 'umap', label = T, label.size = 6)
f613 = FeaturePlot(liver, c('Afp', 'Alb', 'Hnf4a', 'Prox1'))
f1214 = FeaturePlot(liver, c('Cd68', 'Marco', 'Ppbp', 'Itga2b'))
fe = FeaturePlot(liver, c('Hba-a2', 'Hba-a1', 'Hbb-bs', 'Hbb-bt'))
fep = FeaturePlot(liver, c('Klf1', 'Alad', 'Blvrb', 'Gata1'))
ggsave('d.png', d)
ggsave('f613.png', f613, scale = 2)
ggsave('f1214.png', f1214, scale = 2)
ggsave('fe.png', fe, scale = 2)
ggsave('fep.png', fep, scale = 2)

liver10 = subset(liver, subset = seurat_clusters == 10)
if (!file.exists('liver10.rds')){
   liver10 = matrix_processing(liver10, 0.6)
   saveRDS(liver10, 'liver10.rds', compress = T)
} else {
    liver10 = readRDS('liver10.rds')
}
d10 = DimPlot(liver10, reduction = 'umap', label = T, label.size = 6)
cell_Neutr = names(liver10@active.ident[liver10$seurat_clusters == 2])
cell_Stem = names(liver10@active.ident[liver10$seurat_clusters == 0])
f10 = FeaturePlot(liver10, c('Cd34', 'Cmtm7', 'S100a9', 'S100a8'))
f10t = FeaturePlot(liver, c('Cd34', 'Cmtm7', 'S100a9', 'S100a8'))
ggsave('d10.png', d10)
ggsave('f10.png', f10, scale = 2)
ggsave('f10t.png', f10t, scale = 2)


liver11 = subset(liver, subset = seurat_clusters == 11)
if (!file.exists('liver11.rds')){
   liver11 = matrix_processing(liver11, 0.6)
   saveRDS(liver11, 'liver11.rds', compress = T)
} else {
    liver11 = readRDS('liver11.rds')
}
d11 = DimPlot(liver11, reduction = 'umap', label = T, label.size = 6)
cell_Endo = names(liver11@active.ident[liver11$seurat_clusters == 0])
cell_Mesen = names(liver11@active.ident[liver11$seurat_clusters == 1])
f11 = FeaturePlot(liver11, c('Lyve1', 'Kdr', 'Pdgfra', 'Col1a2'))
f11t = FeaturePlot(liver, c('Lyve1', 'Kdr', 'Pdgfra', 'Col1a2'))
ggsave('d11.png', d11)
ggsave('f11.png', f11, scale = 2)
ggsave('f11t.png', f11t, scale = 2)

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

back.up = liver@active.ident
ident = as.character(back.up)
ident = ifelse(ident %in% names(map), map[ident], ident)
ident = ifelse(names(back.up) %in% cell_Endo, 'Endothelial', ident)
ident = ifelse(names(back.up) %in% cell_Mesen, 'Mesenchymal', ident)
ident = ifelse(names(back.up) %in% cell_Neutr, 'Neutrophil', ident)
ident = ifelse(names(back.up) %in% cell_Stem, 'Hematopoietic', ident)
ident = as.factor(ident)
names(ident) = names(back.up)
liver@active.ident = ident
liver$seurat_clusters = ident
saveRDS(liver, 'liver_labeled.rds', compress = T)
dt = DimPlot(liver, reduction = 'umap', label = T)
ggsave('dt.png', dt)

