import matplotlib.pyplot as plt
import sys

f = open(sys.argv[1])
n = [100 * (2 ** k) for k in range(11)]
t = []

for i in range(11):
	t.append(float(f.readline().strip()))

fig = plt.figure(figsize=(8,4))
plt.plot(n, t)
plt.xlabel('n')
plt.ylabel('time (s)')
plt.savefig('qiita/time.png');

fig = plt.figure(figsize=(8,4))
plt.plot(list(range(11)), t)
plt.xlabel('k (n = 100 * 2^k)')
plt.ylabel('time (s)')
plt.savefig('qiita/time-pow.png');
