# suppressPackageStartupMessages(library(gmodels))
x = c(rnorm(5), rnorm(5) + 4)
y = 3 * c(rnorm(5), rnorm(5) + 4)
d = rbind(x, y, a = 0.1 * x, b = 0.2 * x, c = 0.3 * x, o = 0.1 * y, p = 0.2 * y, q = 0.3 * y)
# colnames(d) = paste('s', 1:10, sep = "")
# pca = fast.prcomp(t(dat))
# summary(pca)$importance
# biplot(pca, cex = c(1.3, 1.2))

deg = d
pca = fast.prcomp(t(dat))

