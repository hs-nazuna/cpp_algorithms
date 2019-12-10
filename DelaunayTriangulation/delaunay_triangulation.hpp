#ifndef TSUKASA_DIARY_DELAUNAY_TRIANGULATION_HPP
#define TSUKASA_DIARY_DELAUNAY_TRIANGULATION_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include <limits>
#include <cassert>

namespace delaunay {

static const double PI = acos(-1);

struct Point2D {
	double x, y;
};

struct Edge2D {
	Point2D a, b;
};

struct Triangle2D {
	Point2D a, b, c;
};

Triangle2D get_bounding_triangle(const std::vector< Point2D >& P) {	
	double min_x = std::numeric_limits< double >::max();
	double min_y = std::numeric_limits< double >::max();
	double max_x = std::numeric_limits< double >::min();
	double max_y = std::numeric_limits< double >::min();
	
	size_t n = P.size();
	for (size_t i=0; i<n; ++i) {
		const Point2D& p = P[i];
		min_x = std::min(min_x, p.x);
		min_y = std::min(min_y, p.y);
		max_x = std::max(max_x, p.x);
		max_y = std::max(max_y, p.y);
	}
	
	double cx = (max_x + min_x) / 2;
	double cy = (max_y + min_y) / 2;
	double r = std::max(max_x - min_x, max_y - min_y);
	double dx = 2 * r * cos(PI * 30 / 180);
	double dy = r;
	
	Point2D a{cx, cy + 2 * r}, b{cx - dx, cy - dy}, c{cx  + dx, cy - dy};
	return Triangle2D{a, b, c};
}

std::vector< Edge2D > delaunay_core(const std::vector< Point2D >& P) {
	Triangle2D boundary = get_bounding_triangle(P);
	
	size_t n = P.size();
	std::vector< size_t > id(n);
	
	std::iota(id.begin(), id.end(), 0);
	std::random_device seed_gen;
	std::mt19937 engine(seed_gen());
	std::shuffle(id.begin(), id.end(), engine);
	
	for (size_t i : id) {
		const Point2D& p = P[i];
	}
	
	std::vector< Edge2D > ret;
	return ret;
}

std::vector< Edge2D > delaunay_triangulation(const std::vector< Point2D >& P) { return delaunay_core(P); }
std::vector< Edge2D > delaunay_triangulation(const std::vector< double >& x, const std::vector< double >& y) {
	size_t n = x.size();
	assert(n == y.size());
	std::vector< Point2D > P(n);
	for (size_t i=0; i<n; ++i) {
		P[i].x = x[i];
		P[i].y = y[i];	
	}
	return delaunay_core(P);
}

}

#endif
