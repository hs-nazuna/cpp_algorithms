import sys
import math

PI = math.pi
n, k = int(sys.argv[1]), int(sys.argv[2])

print(k * n + 1)
print(0, 0)

for a in range(k):
	r = a + 1
	for i in range(n):
		print(math.cos(2 * PI * i / n) * r, math.sin(2 * PI * i / n) * r)
