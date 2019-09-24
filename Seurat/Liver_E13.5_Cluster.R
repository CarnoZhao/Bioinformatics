library(dplyr)
library(Seurat)
library(ggplot2)

path = '/mnt/d/Codes/DataLists/GeneMatries/LiverMatrix/'
data = Read10X(data.dir = path)
liver = CreateSeuratObject(counts = data, project = 'Liver', min.cells = 3, min.features = 200)
liver[["percent.mt"]] = PercentageFeatureSet(liver, pattern = "^MT-")
liver = subset(liver, subset = nFeature_RNA > 500 & percent.mt < 10)

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
# 
### first clustering
liver = matrix_processing(liver, 0.8)



liver.markers = FindAllMarkers(liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top12 = liver.markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_logFC)
graph = DoHeatmap(liver, features = top12$gene) + NoLegend()
plot(graph)

DimPlot(liver, reduction = 'umap', label = T, label.size = 5)
DimPlot(liver, reduction = 'tsne', label = T, label.size = 5)

backup_cluster = liver@active.ident
map = c('12' = 'Endothe', '9' = 'Macroph', '10' = 'Hematopo', '4' = 'Hepatocyte', '8' = 'Hepatoblast', '11' = 'Mesenchy', '13' = 'Megakaryocyte')

ident = as.character(liver@active.ident)
ident = ifelse(ident %in% names(map), map[ident], ident)
ident = as.factor(ident)
names(ident) = names(backup_cluster)
liver@active.ident = ident

subliver = subset(liver, subset = seurat_clusters %in% c(0:3, 5:7))
subliver = matrix_processing(subliver, 0.25)
DimPlot(subliver, reduction = 'umap', label = T)
submap = c('0' = 'No signif marker', '1' = 'Eryth_Hba-a2_high', '3' = 'Eryth_Hba-a2_high', '2' = 'Eryth_Klf_high', '4' = 'Mast cell', '5' = 'Mast cell', '6' = 'Mast cell')
subback_up = subliver@active.ident
ident = as.character(subliver@active.ident)
ident = ifelse(ident %in% names(submap), submap[ident], ident)
ident = as.factor(ident)
names(ident) = names(subback_up)
subliver@active.ident = ident
DimPlot(subliver, reduction = 'umap', label = T)
ggsave('/mnt/d/Codes/DataLists/GeneMatries/Layer2.png')

ident = as.character(liver@active.ident)
cells = names(liver@active.ident)
subcells = names(subliver@active.ident)
map = as.character(subliver@active.ident)
names(map) = subcells
ident = ifelse(cells %in% subcells, map[cells], ident)
ident = as.factor(ident)
names(ident) = cells
liver@active.ident = ident
DimPlot(liver, reduction = 'umap', label = T, label.size = 3)
ggsave('/mnt/d/Codes/DataLists/GeneMatries/TotalCluster.png')

liver = readRDS('/mnt/d/Codes/DataLists/GeneMatries/Liver.rds')
markers = FindAllMarkers(liver, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
markers = read.csv('/mnt/d/Codes/DataLists/GeneMatries/markers.csv', row.names = 1, header = T)
top12 = markers %>% group_by(cluster) %>% top_n(n = 12, wt = avg_logFC)
DoHeatmap(liver, features = top12$gene) + NoLegend()


FeaturePlot(liver, c('Lyve1', 'Kdr')) # 12
ggsave('/mnt/d/Codes/DataLists/GeneMatries/marker12.png')

FeaturePlot(liver, c('C1qa', 'C1qb', 'C1qc', 'Cd68')) # 9
ggsave('/mnt/d/Codes/DataLists/GeneMatries/marker9.png')

FeaturePlot(liver, c('Cd34')) # 10
ggsave('/mnt/d/Codes/DataLists/GeneMatries/marker10.png')

FeaturePlot(liver, c('Afp', 'Alb', 'Hnf4a', 'Prox1')) # 48
ggsave('/mnt/d/Codes/DataLists/GeneMatries/marker48.png')

FeaturePlot(liver, c('Vim', 'Mest', 'Pdgfa', 'Col1a2')) # 11
ggsave('/mnt/d/Codes/DataLists/GeneMatries/marker11.png')

FeaturePlot(liver, c('Gp9', 'Myl9', 'Ppbp')) # 13
ggsave('/mnt/d/Codes/DataLists/GeneMatries/marker13.png')

### Quality check

par(mfrow = c(2, 2))
percent = apply(as.matrix(clusters), 1,function(x) {
  subcells = names(liver@active.ident[liver@active.ident == x])
  cnts = subset(liver, cells = subcells)$nCount_RNA
  hist(cnts, xlab = x, main = NULL)
  hist(cnts[cnts < 4000], breaks = 10, xlab = x, main = NULL)
  print(paste(x, ':', sep = ''))
  print(length(cnts[cnts < 2000]) / length(cnts))
})

cnts = data.frame(cnt = liver$nCount_RNA, clusters = liver@active.ident)
ggplot(cnts, aes(y = log10(cnt), x = clusters, color = clusters, fill = clusters)) + geom_boxplot()
cnts$logcnt = log10(cnts$cnt)
t.test(x = cnts$logcnt[cnts$clusters == 'Hepatoblast'],
       y = cnts$logcnt[cnts$clusters == 'No significant marker'])
