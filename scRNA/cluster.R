library(igraph)
library(RANN)
library(ggplot2)
library(easyGgplot2)


d = read.csv('D:/Codes/DataLists/Output/PC1_50.csv', header = T, row.names = 1)
knn = nn2(d, k = 30)
idx = knn$nn.idx
egs = cbind(rep(1:nrow(idx), ncol(idx)), as.vector(idx))
graph = graph_from_edgelist(egs, directed = FALSE)
clst = cluster_louvain(graph)
member = clst$membership

d$member = as.factor(member)
write.csv(member, row.names = F, file = 'D:/Codes/DataLists/Output/member.csv')

g.list = lapply(1:9, function(i){
    g = ggplot(d, aes(x = d[,2 * i - 1], y = d[,2 * i], color = member)) + 
        geom_point(size = 0.3) + 
        theme(legend.position = "none") +
        xlab(paste('PC', 2 * i - 1, sep = '')) + 
        ylab(paste('PC', 2 * i, sep = ''))
    return(g)
})
do.call('ggplot2.multiplot', c(g.list, cols = 3))
