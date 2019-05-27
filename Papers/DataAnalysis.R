library(ggplot2)

d = read.csv('/mnt/d/Codes/DataLists/lncRNA/len_gc.csv')
colnames(d) = c('name', 'length', 'GC')

d = d[!(d$name != 'mRNA' & d$length < 200),]

d$log_len = log10(d$length)
hist1 = ggplot(d[d$name == 'mRNA', ], aes(x = log_len)) + geom_histogram(bins = 50)
plot(hist1)
hist2 = ggplot(d[d$name != 'mRNA', ], aes(x = log_len)) + geom_histogram(bins = 50)
plot(hist2)

d$is_mRNA = d$name == 'mRNA'

print(summary(d[d$name == 'mRNA', 2]))
print(summary(d[d$name != 'mRNA', 2]))

# graph = ggplot(d, aes(x = GC, y = log_len, color = is_mRNA)) + 
#         geom_point() + 
#         xlab('GC%') + 
#         ylab('length')
# plot(graph)

head(d)
