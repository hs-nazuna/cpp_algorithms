#ifndef TSUKASA_DIARY_DELAUNAY_TRIANGULATION_HPP
#define TSUKASA_DIARY_DELAUNAY_TRIANGULATION_HPP

#include <iostream>
#include <vector>
#include <list>
#include <set>
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

class DelaunayTriangulation {	
	/*** Geometry ***/
	struct Point {
		/*** Point on the 2D-plane ***/
		double x, y;
		Point operator + (const Point& p) const { return Point{x + p.x, y + p.y}; }
		Point operator - (const Point& p) const { return Point{x - p.x, y - p.y}; }
		Point operator * (double d) const { return Point{x * d, y * d}; }
		Point operator / (double d) const { return Point{x / d, y / d}; }
		Point normal() const { return Point{y, -x}; }
		double norm() const { return x * x + y + y; }
	};
	
	double cross(Point a, Point b) const { return a.x * b.y - a.y * b.x; }
	bool is_ccw(size_t a, size_t b, size_t c) const { return cross(P[b] - P[a], P[c] - P[a]) > 0; }
	
	using Edge = std::pair<size_t, size_t>;
	Edge make_edge(size_t a, size_t b) const { return Edge(std::min(a,b), std::max(a,b)); }
	
	struct Triangle {
		/*** Triangle on the 2D-plane ***/
		size_t a, b, c;
		
		size_t opposite(Edge e) const {
			if (e.first != a and e.second != a) return a;
			if (e.first != b and e.second != b) return b;
			return c;
		}
	};
	
	Triangle make_triangle(size_t a, size_t b, size_t c) const {
		/*** Make triangle where the direction 'a->b->c' is counter-clockwise ***/
		if (!is_ccw(a, b, c)) std::swap(b, c);
		return Triangle{a, b, c};
	}
	
	bool point_in_triangle(size_t tar, const Triangle& t) const {
		if (!is_ccw(t.a, t.b, tar)) return false;
		if (!is_ccw(t.b, t.c, tar)) return false;
		if (!is_ccw(t.c, t.a, tar)) return false;
		return true;
	}
	
private:
	size_t n;
	std::vector<Point> P;
	std::vector<Triangle> T;
	std::vector<std::list<size_t>> CH;
	
	size_t rdd_add_child(size_t k, const Triangle& t) {
		size_t new_k = T.size();
		T.push_back(t);
		CH.push_back(std::list<size_t>());
		CH[k].push_back(new_k);
		return new_k;
	}
	
	size_t rdd_add_child(size_t k1, size_t k2, const Triangle& t) {
		size_t new_k = T.size();
		T.push_back(t);
		CH.push_back(std::list<size_t>());
		CH[k1].push_back(new_k);
		CH[k2].push_back(new_k);
		return new_k;
	}
	
	size_t find_child(size_t k, size_t tar) const {
		for (size_t i : CH[k]) if (point_in_triangle(tar, T[i])) return i;
		return std::numeric_limits<size_t>::max();
	}
	
private:	
	struct EdgeHash {
		size_t operator () (const Edge& e) const {
			static const size_t FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();
			size_t key = (e.first << 16) ^ e.second;
			key ^= FIXED_RANDOM;
			return key ^ (key >> 16);
		}
	};
	
	std::unordered_map<Edge, std::set<size_t>, EdgeHash> e2t;
	
	void register_triangle(size_t k, const Triangle& t) {
//		std::cerr << "Register " << t.a << " " << t.b << " " << t.c << std::endl;
		e2t[make_edge(t.a, t.b)].insert(k);
		e2t[make_edge(t.b, t.c)].insert(k);
		e2t[make_edge(t.c, t.a)].insert(k);
	}
	
	void unregister_triangle(size_t k, const Triangle& t) {
//		std::cerr << "Unregister " << t.a << " " << t.b << " " << t.c << std::endl;
		e2t[make_edge(t.a, t.b)].erase(k);
		e2t[make_edge(t.b, t.c)].erase(k);
		e2t[make_edge(t.c, t.a)].erase(k);		
	}
	
private:
	bool is_ilegal(Edge e, size_t c, size_t tar) {
		size_t a = e.first, b = e.second;
		
		Point s1 = (P[a] + P[b]) / 2.;
		Point t1 = s1 + (P[b] - P[a]).normal();
		
		Point s2 = (P[b] + P[c]) / 2.;
		Point t2 = s2 + (P[c] - P[b]).normal();
		
		double d1 = cross(t2 - s2, s2 - s1);
		double d2 = cross(t2 - s2, t1 - s1);
		Point center = s1 + (t1 - s1) * d1 / d2;
		
		double x = center.x, y = center.y;
			
		{
			double dx1 = P[a].x - x, dy1 = P[a].y - y;
			double dx2 = P[b].x - x, dy2 = P[b].y - y;
			double dx3 = P[c].x - x, dy3 = P[c].y - y;
//			std::cerr << (dx1 * dx1 + dy1 * dy1) << " " << (dx2 * dx2 + dy2 * dy2) << std::endl;
//			std::cerr << (dx2 * dx2 + dy2 * dy2) << " " << (dx3 * dx3 + dy3 * dy3) << std::endl;
//			std::cerr << (dx3 * dx3 + dy3 * dy3) << " " << (dx1 * dx1 + dy1 * dy1) << std::endl;
			assert(std::abs((dx1 * dx1 + dy1 * dy1) - (dx2 * dx2 + dy2 * dy2)) < 1e-5);
			assert(std::abs((dx2 * dx2 + dy2 * dy2) - (dx3 * dx3 + dy3 * dy3)) < 1e-5);
			assert(std::abs((dx3 * dx3 + dy3 * dy3) - (dx1 * dx1 + dy1 * dy1)) < 1e-5);
		}
		
		double dx1 = P[a].x - x, dy1 = P[a].y - y;
		double dx2 = P[tar].x - x, dy2 = P[tar].y - y;
		return (dx1 * dx1 + dy1 * dy1) > (dx2 * dx2 + dy2 * dy2);
	}

	void legalize_edge(size_t piv, Edge e) {
		if (e2t.count(e) == 0) return;
		if (e2t[e].size() != 2) return;
		
		size_t k1 = *e2t[e].begin();
		size_t k2 = *e2t[e].rbegin();
		size_t a = T[k1].opposite(e);
		size_t b = T[k2].opposite(e);
		
		if (is_ilegal(e, a, b)) {
			unregister_triangle(k1, T[k1]);
			unregister_triangle(k2, T[k2]);
			e2t.erase(e);
			
			Triangle t1 = make_triangle(e.first, a, b);
			Triangle t2 = make_triangle(e.second, a, b);
			
			size_t k1_ = rdd_add_child(k1, k2, t1);
			size_t k2_ = rdd_add_child(k1, k2, t2);
			
			register_triangle(k1_, t1);
			register_triangle(k2_, t2);
			
			size_t c = (piv != a ? a : b);
			legalize_edge(piv, make_edge(e.first, c));
			legalize_edge(piv, make_edge(e.second, c));
		}
	}
	
	void interior_division(size_t k, size_t tar) {
		Triangle t1 = make_triangle(T[k].a, T[k].b, tar);
		Triangle t2 = make_triangle(T[k].b, T[k].c, tar);
		Triangle t3 = make_triangle(T[k].c, T[k].a, tar);
		
		unregister_triangle(k, T[k]);
		
		size_t k1 = rdd_add_child(k, t1);
		size_t k2 = rdd_add_child(k, t2);
		size_t k3 = rdd_add_child(k, t3);
		
		register_triangle(k1, t1);
		register_triangle(k2, t2);
		register_triangle(k3, t3);
		
		legalize_edge(tar, make_edge(T[k].a, T[k].b));
		legalize_edge(tar, make_edge(T[k].b, T[k].c));
		legalize_edge(tar, make_edge(T[k].c, T[k].a));
	}
	
	void edge_division(size_t k, size_t tar) {
		std::cerr << "Edge division @" << k << std::endl;
	}
	
public:
	DelaunayTriangulation(const std::vector<double>& x, const std::vector<double>& y) {
		n = x.size();
		P = std::vector<Point>(n);
		for (size_t i=0; i<n; ++i) {
			P[i].x = x[i];
			P[i].y = y[i];
		}
		
		double r = std::numeric_limits<double>::min();
		for (size_t i=0; i<n; ++i) {
			r = std::max(r, std::abs(P[i].x));
			r = std::max(r, std::abs(P[i].y));
		}
		
		P.push_back(Point{3.1 * r, 0});
		P.push_back(Point{0, 3.1 * r});
		P.push_back(Point{-3.1 * r, -3.1 * r});
	}
	
	void core() {
		/*** Core algorithm ***/
		std::vector<size_t> id(n);
		std::iota(id.begin(), id.end(), size_t(0));
		std::random_device seed_gen;
		std::mt19937 engine(seed_gen());
		std::shuffle(id.begin(), id.end(), engine);
		
		T = std::vector<Triangle>();
		CH = std::vector<std::list<size_t>>();
		T.push_back(make_triangle(n, n+1, n+2));
		CH.push_back(std::list<size_t>());
		register_triangle(0, T[0]);
		e2t = std::unordered_map<Edge, std::set<size_t>, EdgeHash>();
		
		for (size_t tar : id) {
			size_t k = 0;
			while (!CH[k].empty()) {
				size_t nxt = find_child(k, tar);
				if (nxt == std::numeric_limits<size_t>::max()) break;
				k = nxt;
			}
			if (CH[k].empty()) interior_division(k, tar);
			else edge_division(k, tar);
		}
	}
	
	void dump(std::ostream& os) {
		os << P[n].x << " " << P[n].y << std::endl;
		os << P[n+1].x << " " << P[n+1].y << std::endl;
		os << P[n+2].x << " " << P[n+2].y << std::endl;
		os << e2t.size() << std::endl;
//		std::cerr << "#Edge: " << e2t.size() << std::endl;
		for (auto it : e2t) {
			Edge e = it.first;
			os << e.first << " " << e.second << std::endl;
//			std::cerr << e.first << "-" << e.second << " x " << it.second.size() << std::endl;
		}
	}
};

}

#endif
