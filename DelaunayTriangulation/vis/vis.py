import sys
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import numpy as np

# Read input file
f = open(sys.argv[1])
n = int(f.readline())
x, y = list(), list()
for i in range(n):
    xi, yi = map(float, f.readline().split())
    x.append(xi)
    y.append(yi)

# Read output file
f = open(sys.argv[2])
m = int(f.readline())
edge = list()
for i in range(m):
    a, b = map(int, f.readline().split())
    edge.append([(x[a], y[a]), (x[b], y[b])])

# Visualize
fig = plt.figure(figsize=(4,4))
ax = fig.add_axes([0.1,0.1,0.85,0.85])#add_subplot(1,1,1)
xy_max = max(np.max(x), np.max(y)) + 1
xy_min = min(np.min(x), np.min(y)) - 1
ax.set_xlim(xy_min, xy_max)
ax.set_ylim(xy_min, xy_max)
ax.scatter(x, y, s=20)
ax.add_collection(mc.LineCollection(edge))
plt.savefig(sys.argv[3])
