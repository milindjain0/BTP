import numpy as np
from sklearn.cluster import KMeans
file1 = open("l2.txt")
data = []
count = 0
for line in file1:
    if count == -1:
        count+=1
        continue
    line = line.replace("\n","")
    line = line.split(" ")
    data.append([float(i) for i in line])
X = np.array(data)
kmeans = KMeans(n_clusters=2	, random_state=0).fit(X)
print kmeans.cluster_centers_
print kmeans.score(X)
