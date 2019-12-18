test: test.cpp ../delaunay_triangulation.hpp
	g++ -std=c++17 -O3 -o test test.cpp

verify: mst_verify.cpp ../delaunay_triangulation.hpp
	g++ -std=c++17 -O3 -o verify mst_verify.cpp

measure: measure.cpp ../delaunay_triangulation.hpp
	g++ -std=c++17 -O3 -o measure measure.cpp

clean:
	rm test verify measure
