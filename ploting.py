#from matplotlitb import pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

N = 50
x1 = np.random.random_integers(10,70,80)
y1 = np.random.random_integers(20,70,80)

x2 = np.random.random_integers(100,160,200)
y2 = np.random.random_integers(10,150,200)

x3 = np.random.random_integers(90,100,40)
y3 = np.random.random_integers(,100,40)

x = np.append(x1,x2)
y = np.append(y1,y2)
x = np.append(x,x3)
y = np.append(y,y3)

xl = list(x)
yl = list(y)

import sys
for i in range(len(xl)):
    # if i == len(xl)-1:
    #     sys.stdout.write(str(xl[i])+" "+str(yl[i]))
    # else:
    sys.stdout.write(str(xl[i])+" "+str(yl[i])+"\n")


data = []
count = 0
for i in range(len(xl)):
    data.append([xl[i],yl[i]])
X = np.array(data)

kmeans = KMeans(n_clusters=2	, random_state=0).fit(X)
centroids = kmeans.cluster_centers_
print kmeans.cluster_centers_
print kmeans.score(X)
plt.scatter(centroids[:, 0], centroids[:, 1],
            marker='x', s=169, linewidths=3,
            color='yellow', zorder=10)
plt.scatter(x1, y1,color="blue")
plt.scatter(x2, y2,color="red")

plt.show()
