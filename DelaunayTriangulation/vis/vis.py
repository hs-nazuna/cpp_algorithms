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
for i in range(3):
	xi, yi = map(float, f.readline().split())
	x.append(xi)
	y.append(yi)
m = int(f.readline())
edge = list()
for i in range(m):
	a, b = map(int, f.readline().split())
	edge.append([(x[a], y[a]), (x[b], y[b])])

# Visualize
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1,1,1)
cx, cy = np.mean(x), np.mean(y)
delta_x, delta_y = np.max(x) - np.min(x), np.max(y) - np.min(y)
delta = max(delta_x, delta_y)
ax.set_xlim(cx - delta, cx + delta)
ax.set_ylim(cy - delta, cy + delta)
ax.scatter(x, y)
ax.add_collection(mc.LineCollection(edge))
plt.savefig(sys.argv[3])
