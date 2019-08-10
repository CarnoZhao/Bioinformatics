library(Seurat)
library(dplyr)
library(ggplot2)
peak.matrix = Read10X_h5('../../data/scatac/filtered_peak_bc_matrix.h5')
activity.matrix = CreateGeneActivityMatrix(
    peak.matrix = peak.matrix,
    annotation.file = '../../reference/mouse/Mus_musculus.GRCm38.97.gtf.gz',
    seq.levels = c(1:19, 'X', 'Y'),
    upstream = 2000,
    verbose = T
)

liver.atac = CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
liver.atac[["ACTIVITY"]] = CreateAssayObject(counts = activity.matrix)
meta = read.table(
    "../../data/scatac/singlecell.csv", 
    sep = ",", 
    header = T, 
    row.names = 1, 
    stringsAsFactors = F
)
meta = meta[colnames(liver.atac), ]
liver.atac = AddMetaData(liver.atac, metadata = meta)
liver.atac = subset(liver.atac, subset = nCount_ATAC > 5000)
liver.atac$tech = "atac"

DefaultAssay(liver.atac) <- "ACTIVITY"
liver.atac <- FindVariableFeatures(liver.atac)
liver.atac <- NormalizeData(liver.atac)
liver.atac <- ScaleData(liver.atac)
DefaultAssay(liver.atac) <- "ATAC"
VariableFeatures(liver.atac) <- names(which(Matrix::rowSums(liver.atac) > 100))
liver.atac <- RunLSI(liver.atac, n = 50, scale.max = NULL)
liver.atac <- RunUMAP(liver.atac, reduction = "lsi", dims = 1:50)
saveRDS(liver.atac, '../../data/scatac/liver.atac.rds')

liver.atac = readRDS('../Download/liver.atac.rds')
