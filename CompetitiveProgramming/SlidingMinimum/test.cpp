#include <iostream>
#include <vector>
using namespace std;

#include "sliding_minimum.hpp"

int main() {
	int n, K;
	cin >> n >> K;
	
	vector<int> v(n);
	for (int i=0; i<n; ++i) cin >> v[i];
	
	for (int i=0; i<n; ++i) cout << v[i] << (i<n-1 ? " " : "\n");
	cout << "Window Size : " << K << endl;
	
	vector<size_t> sm = sliding_minimum(v, K);
	
	cout << "Index : ";
	for (int i=0; i<n-K+1; ++i) cout << sm[i] << (i<n-K ? " " : "\n");
	cout << "Value : ";
	for (int i=0; i<n-K+1; ++i) cout << v[sm[i]] << (i<n-K ? " " : "\n");
}
