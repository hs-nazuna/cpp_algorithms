#include "../delaunay_triangulation.hpp"

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

int main() {
    int n;
    std::cin >> n;
    
    std::vector<double> x(n), y(n);
    for (int i=0; i<n; ++i) std::cin >> x[i] >> y[i];
    
    Timer timer;
    delaunay::DelaunayTriangulation DT(x, y);
    DT.execute();
    std::cout << timer.curTime() << std::endl;
}
