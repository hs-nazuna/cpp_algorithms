import sys

h, w = int(sys.argv[1]), int(sys.argv[2])
print(h * w)
for y in range(h):
	for x in range(w):
		print(y, x)
