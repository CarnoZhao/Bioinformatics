import tsne
import numpy as np

d = np.genfromtxt('D:/Codes/DataLists/Output/PC1-50.csv', delimiter = ',', skip_header = 1, usecols = range(1, 51))
print(d.shape)

label = np.genfromtxt('D:/Codes/DataLists/Output/member.csv', delimiter = ',', skip_header = 1)
Y = tsne.tsne(d, 2, 50, 30.0)
np.savetxt('D:/Codes/DataLists/Output/tsne/csv', Y, delimiter = ',')
