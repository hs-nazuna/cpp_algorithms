#include <iostream>
using namespace std;

#include "RandomizedTopologicalSort.hpp"

using RTS = RandomizedTopologicalSort;

int main() {
	int n, m;
	cin >> n >> m;
	
	vector< vector< int > > adj_list(n);
	for (int i=0; i<m; ++i) {
		int f, t;
		cin >> f >> t;
		adj_list[f].push_back(t);
	}
	
	RTS rts;
	
	for (int i=0; i<10; ++i) {
		vector< int > order = rts.generate(adj_list);
		for (int v : order) cout << v << " ";
		cout << endl;
	}
}
