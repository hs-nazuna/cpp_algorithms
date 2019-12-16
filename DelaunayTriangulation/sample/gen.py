import random
import sys

n = int(sys.argv[1])
P = set()
while n > 0:
	p = (round(random.random(), 3), round(random.random(), 3))
	if p in P:
		continue;
	P.add(p)
	n -= 1

print(len(P))
for p in P:
	print(p[0], p[1])
