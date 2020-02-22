import random
import sys

n = int(sys.argv[1])
c_min = int(sys.argv[2])
c_max = int(sys.argv[3])
P = set()
while n > 0:
    x = c_min + (c_max - c_min) * random.random()
    y = c_min + (c_max - c_min) * random.random()
    p = (round(x, 3), round(y, 3))
    if p in P:
        continue;
    P.add(p)
    n -= 1

print(len(P))
for p in P:
    print(p[0], p[1])
