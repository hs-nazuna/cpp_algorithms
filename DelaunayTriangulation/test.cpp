#include <iostream>

#include "delaunay_triangulation.hpp"

int main() {
	int n;
	std::cin >> n;
	
	std::vector< double > x(n), y(n);
	for (int i=0; i<n; ++i) std::cin >> x[i] >> y[i];
	
	delaunay::delaunay_triangulation(x, y);
}
