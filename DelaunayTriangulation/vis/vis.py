import sys
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import numpy as np

# Open File
f = open(sys.argv[1])

# Read Points
n = int(f.readline())
x, y = list(), list()
for i in range(n):
	xi, yi = map(float, f.readline().split())
	x.append(xi)
	y.append(yi)

# Read Edges
m = int(f.readline())
edge = list()
for i in range(m):
	a, b = map(int, f.readline().split())
	edge.append([(x[a], y[a]), (x[b], y[b])])

# Visualize
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
cx, cy = np.mean(x), np.mean(y)
delta = (np.max(np.concatenate([x,y])) - np.min(np.concatenate([x,y]))) / 2
ax.set_xlim(cx - delta, cx + delta)
ax.set_ylim(cy - delta, cy + delta)
ax.scatter(x, y)
ax.add_collection(mc.LineCollection(edge))
plt.savefig(sys.argv[2])
