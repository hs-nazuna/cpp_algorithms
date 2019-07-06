#ifndef UNIFORM_RANDOM_SELECTOR_HPP
#define UNIFORM_RANDOM_SELECTOR_HPP

#include <vector>
#include <cassert>
#include "XorShift.hpp"

template< typename T >
class UniformRandomSelector {
private:
	XorShift rnd;
	
	int max_n;
	int n;
	std::vector< T > data;
	
	void expand() {
		if (n == max_n) {
			max_n *= 2;
			data.resize(max_n);
		}
	}
	
public:
	UniformRandomSelector() : max_n(1), n(0), data(1) {}
	UniformRandomSelector(int n_) : max_n(n_), n(0), data(n_) {}
		
	void add(const T& t) {
		expand();
		data[n++] = t;
	}
	
	T select() {
		assert(n > 0);
		int i = rnd.nextInt(n);
		swap(data[i], data[--n]);
		return data[n];
	}
	
	int size() { return n; }
	
	void clean() { n = 0; }
};

template< typename T > using URS = UniformRandomSelector< T >;

#endif
