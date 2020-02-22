#include "binary_trie.hpp"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <cassert>
#include <sys/time.h>

struct Timer {
	struct timeval start, cur;
	double limit;
	
	Timer() : limit(0) { gettimeofday(&start, NULL); }
	Timer(double l) : limit(l) { gettimeofday(&start, NULL); }
	
	bool isLimit() { return curTime() > limit; }
	
	double curTime() {
		gettimeofday(&cur, NULL);
		return (cur.tv_sec - start.tv_sec) + (cur.tv_usec - start.tv_usec) / 1e6;
	}
};

struct XorShift {
	unsigned int x, y, z, w;
	
	XorShift() { init(); }
	XorShift(unsigned int seed) { init(seed); }
	
	void init(unsigned int seed = 88675123) { x = 123456789; y = 362436069; z = 521288629; w = seed; }
	
	void setNext() {
		unsigned int t = x ^ (x << 11);
		x = y; y = z; z = w;
		w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
	}
	
	unsigned int nextInt() { setNext(); return w; }
	unsigned int nextInt(unsigned int n) { return nextInt() % n; }
	unsigned int nextInt(unsigned int a, unsigned int b) { assert(a <= b); return a + nextInt(b - a + 1); }
	double nextDouble() { return nextInt() * (1.0 / 0xFFFFFFFFu); }
};

int main() {
	using Pair = std::pair<size_t,size_t>;
	struct Hash { size_t operator () (const Pair& p) const { return (p.first << 16) ^ p.second; } };
	
	
	
	Timer timer;
	XorShift random;
	int n = 1000000;
	std::vector< size_t > a, b, x;
	for (int i=0; i<n; ++i) {
		a.push_back(random.nextInt(100000));
		b.push_back(random.nextInt(100000));
		x.push_back(random.nextInt(100000));
	}
	std::cout << "Instance : " << timer.curTime() << std::endl;
	
	
	
	timer = Timer();
	DoubleBinaryTrie<size_t> trie;
	for (int i=0; i<n; ++i) trie.add(a[i], b[i], x[i]);
	for (int i=0; i<n; ++i) trie.erase(a[i], b[i]);
	std::cout << "Trie : " << timer.curTime() << std::endl;
	
	
	
	timer = Timer();
	std::unordered_map< Pair, std::vector<size_t>, Hash> map;
	for (int i=0; i<n; ++i) map[Pair(a[i],b[i])].push_back(x[i]);
	for (int i=0; i<n; ++i) map.erase(Pair(a[i],b[i]));
	std::cout << "unordered_map : " << timer.curTime() << std::endl;
	
	
	
	timer = Timer();
	std::map< Pair, std::vector<size_t>> map2;	
	for (int i=0; i<n; ++i) map2[Pair(a[i],b[i])].push_back(x[i]);
	for (int i=0; i<n; ++i) map2.erase(Pair(a[i],b[i]));
	std::cout << "map : " << timer.curTime() << std::endl;
}
