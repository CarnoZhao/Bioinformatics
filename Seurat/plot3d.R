library(rgl)
umap.atac = read.csv('D:/Codes/Bioinformatics/Seurat/umap.atac.csv')
umap.rna = read.csv('D:/Codes/Bioinformatics/Seurat/umap.rna.csv')
colnames(umap.rna) = c('type', 'x', 'y', 'z')
colnames(umap.atac) = c('type', 'x', 'y', 'z')

types = names(table(umap.rna$type))
mapping = 1:length(types)
names(mapping) = types

p = function(){
    plot3d(
        c(umap.rna$x, umap.atac$x), 
        c(umap.rna$y, umap.atac$y), 
        c(umap.rna$z, umap.atac$z),
        col = mapping[c(umap.rna$type, umap.atac$type)],
        xlab = "",
        ylab = "",
        zlab = "",
        box = T, axes = T
        )
    
    points.list = lapply(types, function(t) {
        rna = umap.rna[umap.rna$type == t,]
        atac = umap.atac[umap.atac$type == t,]
        num.lines = round(0.05 * min(nrow(rna), nrow(atac)))
        rna = rna[sample(1:nrow(rna), num.lines),]
        atac = atac[sample(1:nrow(rna), num.lines),]
        cbind(rna$x, atac$x, rna$y, atac$y, rna$z, atac$z, t)
    })
    lines = do.call(rbind, points.list)
    apply(lines, MARGIN = 1, function(line){
        xs = line[1:2]
        ys = line[3:4]
        zs = line[5:6]
        type = line[7]
        lines3d(xs, ys, zs, add = T, col = mapping[type], alpha = 0.2)
    })
}

p()

umap.co = read.csv('D:/Codes/Bioinformatics/Seurat/umap.co.csv')
colnames(umap.co) = c('type', 'x', 'y', 'z')
types = names(table(umap.co$type))
mapping = 1:length(types)
names(mapping) = types

p2 = function(){
    plot3d(
        umap.co$x, 
        umap.co$y, 
        umap.co$z,
        col = mapping[umap.co$type],
        xlab = "",
        ylab = "",
        zlab = "",
        box = T, axes = T
    )
    
    points.list = lapply(types, function(t) {
        rna = umap.co[umap.co$type == t & umap.co$z == 1,]
        atac = umap.co[umap.co$type == t & umap.co$z == 0,]
        num.lines = round(0.1 * min(nrow(rna), nrow(atac)))
        rna = rna[sample(1:nrow(rna), num.lines),]
        atac = atac[sample(1:nrow(rna), num.lines),]
        cbind(rna$x, atac$x, rna$y, atac$y, rna$z, atac$z, t)
    })
    lines = do.call(rbind, points.list)
    apply(lines, MARGIN = 1, function(line){
        xs = line[1:2]
        ys = line[3:4]
        zs = line[5:6]
        type = line[7]
        lines3d(xs, ys, zs, add = T, col = mapping[type], alpha = 0.1)
    })
}
p2()
