#ifndef TSUKASA_DIARY_DELAUNAY_TRIANGULATION_HPP
#define TSUKASA_DIARY_DELAUNAY_TRIANGULATION_HPP

#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>
#include <cassert>

namespace delaunay {

static const double PI = acos(-1);

template<typename T> struct Point {
	/*** Point on the 2D-plane ***/
	T x, y;
	size_t id;
	
	T cross(const Point<T>& a) const { return x * a.y - y * a.x; }
	Point<T> operator - (const Point<T>& a) const { return Point<T>{x - a.x, y - a.y}; }
};

template<typename T> bool is_ccw(const Point<T>& a, const Point<T>& b, const Point<T>& c) {
	/*** Check whether the direction 'a->b->c' is counter-clockwise or not ***/
	Point<T> ab = b - a;
	Point<T> ac = c - a;
	return ab.cross(ac) > 0;
}

template<typename T> struct Triangle {
	/*** Triangle with three points (vertices) 'a', 'b', and 'c' on the 2D-plane ***/
	Point<T> a, b, c;
	
	void ccw_arange() {
		// Make the direction 'a->b->c' counter-clockwise
		if (!is_ccw(a, b, c)) std::swap(b, c);
	}
	
	bool has_point(const Point<T>& p) const {
		if (!is_ccw(a, b, p)) return false;
		if (!is_ccw(b, c, p)) return false;
		if (!is_ccw(c, a, p)) return false;
		return true;
	}
};

template<typename T> struct RDD_Node {
	/*** Node of Region Decision Diagrams (RDDs) ***/
	Triangle<T> t;
	std::vector<size_t> ch;
	
	bool is_leaf() const { return ch.empty(); }
};

template<typename T> struct RDD {
	/*** Region Decision Diagram (RDD) ***/
	std::vector<RDD_Node<T>> node;
	
	void initialize(Triangle<T> root_t) {
		RDD_Node<T> root{root_t, std::vector<size_t>()};
		node = std::vector<RDD_Node<T>>(1, root);
	}
	
	size_t add_child(size_t k, const Triangle<T>& t) {
		size_t i = node.size();
		node[k].ch.push_back(i);
		node.push_back(RDD_Node<T>{t, std::vector<size_t>()});
		return i;
	}
	
	size_t find_child(size_t k, const Point<T>& p) const {
		for (size_t i : node[k].ch) if (node[i].t.has_point(p)) return i;
		return std::numeric_limits<size_t>::max();
	}
	
	const RDD_Node<T>& operator [] (int i) { return node[i]; }
};

using Edge = std::pair<size_t, size_t>;

struct EdgeHash {
	size_t operator () (const Edge& e) const {
		static const size_t FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();
		size_t key = (e.first << 16) ^ e.second;
		key ^= FIXED_RANDOM;
		return key ^ (key >> 16);
	}
};

using Edge2Triangles = std::unordered_map<Edge, std::vector<size_t>, EdgeHash>;

Edge make_edge(size_t a, size_t b) {
	if (a > b) std::swap(a, b);
	return Edge(a, b);
}

template<typename T> void register_triangle(Edge2Triangles& e2t, const Triangle<T>& t, size_t k) {
	e2t[make_edge(t.a.id, t.b.id)].push_back(k);
	e2t[make_edge(t.b.id, t.c.id)].push_back(k);
	e2t[make_edge(t.c.id, t.a.id)].push_back(k);
}

template<typename T> void interior_division(Edge2Triangles& e2t, RDD<T>& rdd, size_t k, const Point<T>& p) {
	/*** Subdivide when p lies in the interior of rdd[k].t ***/
	const Triangle<T>& t = rdd[k].t;
	
	Triangle<T> t1{t.a, t.b, p};
	Triangle<T> t2{t.b, t.c, p};
	Triangle<T> t3{t.c, t.a, p};
	
	t1.ccw_arange();
	t2.ccw_arange();
	t3.ccw_arange();
	
	size_t k1 = rdd.add_child(k, t1);
	size_t k2 = rdd.add_child(k, t2);
	size_t k3 = rdd.add_child(k, t3);
	
	register_triangle(e2t, t1, k1);
	register_triangle(e2t, t2, k2);
	register_triangle(e2t, t3, k3);
	
	// TODO : legalize_edge method
}

template<typename T> void edge_division(Edge2Triangles& e2t, RDD<T>& rdd, size_t k, const Point<T>& p) {
	/*** Subdivide when p falls on an edge between two adjacent triangles in rdd[k].t.ch ***/	
	
}

template<typename T> std::vector<Edge> delaunay_core(const std::vector<Point<T>>& P) {
	/*** Main process of the delaunay triangulation ***/	
	size_t n = P.size();
	std::vector<size_t> id(n);
	std::iota(id.begin(), id.end(), 0);
	
	std::random_device seed_gen;
	std::mt19937 engine(seed_gen());
	std::shuffle(id.begin(), id.end(), engine);
	
	Edge2Triangles e2t;
	RDD<T> rdd;
	
	T r = std::numeric_limits<T>::min();
	for (const Point<T>& p : P) {
		r = std::max(r, std::abs(p.x));
		r = std::max(r, std::abs(p.y));
	}
	Point<T> a{4 * r, 0};
	Point<T> b{0, 4 * r};
	Point<T> c{-4 * r, -4 * r};
	rdd.initialize(Triangle<T>{a,b,c});
	
	for (size_t i : id) {
		const Point<T>& p = P[i];		
		size_t k = 0;
		while (!rdd[k].is_leaf()) {
			size_t nxt = rdd.find_child(k, p);
			if (nxt == std::numeric_limits<size_t>::max()) break;
			k = nxt;
		}
		std::cerr << k << std::endl; // debug
		if (rdd[k].is_leaf()) interior_division(e2t, rdd, k, p);
		else edge_division(e2t, rdd, k, p);
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
