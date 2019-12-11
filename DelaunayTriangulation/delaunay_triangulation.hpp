#ifndef TSUKASA_DIARY_DELAUNAY_TRIANGULATION_HPP
#define TSUKASA_DIARY_DELAUNAY_TRIANGULATION_HPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <cmath>
#include <limits>
#include <cassert>

namespace delaunay {

using Edge = std::pair<size_t, size_t>;

static const double PI = acos(-1);

template<typename T> struct Point {
	/* Point on the 2D-plane */
	T x, y;
	size_t id;
	
	T cross(const Point<T>& a) const { return x * a.y - y * a.x; }
	Point<T> operator - (const Point<T>& a) const { return Point<T>{x - a.x, y - a.y}; }
};

template<typename T> bool is_ccw(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
	/* Check whether the direction 'a->b->c' is counter-clockwise or not */
	Point<T> ab = b - a, ac = c - a;
	return ab.cross(ac) > 0;
}

template<typename T> struct Triangle {
	/* Triangle with three points (vertices) 'a', 'b', and 'c' on the 2D-plane */
	Point<T> a, b, c;
	
	void ccw_arrange() { if (!is_ccw(a, b, c)) std::swap(b, c); }
	
	bool has_point(const Point<T>& p) const {
		if (!is_ccw(a, b, p)) return false;
		if (!is_ccw(b, c, p)) return false;
		if (!is_ccw(c, a, p)) return false;
		return true;
	}
};

template<typename T> struct RDD_Node {
	/* Node of Region Decision Diagrams (RDDs) */
	Triangle<T> t;
	std::vector<size_t> ch;
};

template<typename T> struct RDD {
	/* Region Decision Diagram (RDD) */
	std::vector<RDD_Node<T>> node;
	
	void initialize(Triangle<T> root_t) {
		RDD_Node<T> root{root_t, std::vector<size_t>()};
		node = std::vector<RDD_Node<T>>(1, root);
	}
	
	size_t find_leaf(const Point<T>& p) const {
		size_t i = 0;
		while (!node[i].ch.empty()) {
			for (size_t j : node[i].ch) {
				if (node[j].t.has_point(p)) {
					i = j;
					break;
				}
			}
		}
		return i;
	}
};

template<typename T> Triangle<T> get_bounding_triangle(const std::vector<Point<T>>& P) {	
	/* Make a bounding triangle of the given point set 'P' on the 2D-plane */
	T min_x = std::numeric_limits<T>::max(), min_y = std::numeric_limits<T>::max();
	T max_x = std::numeric_limits<T>::min(), max_y = std::numeric_limits<T>::min();
	
	for (const Point<T>& p : P) {
		min_x = std::min(min_x, p.x); min_y = std::min(min_y, p.y);
		max_x = std::max(max_x, p.x); max_y = std::max(max_y, p.y);
	}
	
	T cx = (max_x + min_x) / 2, cy = (max_y + min_y) / 2, r = std::max(max_x - min_x, max_y - min_y);
	T dx = 2 * r * cos(PI * 30 / 180), dy = r;
	return Triangle<T>{{cx, cy + 2 * r, P.size()}, {cx - dx, cy - dy, P.size() + 1}, {cx + dx, cy - dy, P.size() + 2}};
}

template<typename T> std::vector<Edge> delaunay_core(const std::vector<Point<T>>& P) {
	/* Main process of the delaunay triangulation */	
	size_t n = P.size();
	std::vector<size_t> id(n);
	std::iota(id.begin(), id.end(), 0);
	
	std::random_device seed_gen;
	std::mt19937 engine(seed_gen());
	std::shuffle(id.begin(), id.end(), engine);
	
	RDD<T> rdd;
	rdd.initialize(get_bounding_triangle(P));
	
	for (size_t i : id) {
		const Point<T>& p = P[i];
		size_t node_i = rdd.find_leaf(p);
	}
	
	std::vector<Edge> ret;
	return ret;
}

template<typename T> std::vector<Edge> delaunay_triangulation(const std::vector<T>& x, const std::vector<T>& y) {
	size_t n = x.size();
	assert(n == y.size());
	std::vector<Point<T>> P(n);
	for (size_t i=0; i<n; ++i) {
		P[i].x = x[i];
		P[i].y = y[i];
		P[i].id = i;	
	}
	return delaunay_core(P);
}

}

#endif
