library(dplyr)
library(Seurat)

data = Read10X(data.dir = '/p200/liujiang_group/yinyao/Dataset/Seurat/GeneMatries/PBMC/')
pbmc = CreateSeuratObject(counts = data, project = 'pbmc3k', min.cells = 3, min.features = 200)


