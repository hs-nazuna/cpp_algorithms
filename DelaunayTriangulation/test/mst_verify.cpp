#include "../delaunay_triangulation.hpp"

double dst(double x, double y, double xx, double yy) {
    double dx = x - xx, dy = y - yy;
    return std::sqrt(dx * dx + dy * dy);
}

double euclid_mst(const std::vector<double>& x,
                  const std::vector<double>& y,
                  std::vector<delaunay::Edge>& edge)
{
    std::sort(edge.begin(), edge.end(),
              [&](const delaunay::Edge& a, const delaunay::Edge& b) {
                  double A = dst(x[a.first], y[a.first], x[a.second], y[a.second]);
                  double B = dst(x[b.first], y[b.first], x[b.second], y[b.second]);
                  return A < B;
              });
    
    size_t n = x.size();
    std::vector<std::set<size_t>> uf(n);
    std::vector<size_t> root(n);
    for (size_t i=0; i<n; ++i) {
        uf[i].insert(i);
        root[i] = i;
    }	
    
    double ret = 0;
    for (const std::pair<size_t,size_t>& e : edge) {
        size_t r = root[e.first];
        if (uf[r].count(e.second)) continue;
        
        ret += dst(x[e.first], y[e.first], x[e.second], y[e.second]);
        
        size_t a = root[e.first], b = root[e.second];
        if (uf[a].size() < uf[b].size()) std::swap(a, b);
        
        for (size_t v : uf[b]) root[v] = a;
        uf[a].insert(uf[b].begin(), uf[b].end());
    }
    
    return ret;
}

int main() {
    int n;
    std::cin >> n;
    
    std::vector<double> x(n), y(n);
    for (int i=0; i<n; ++i) std::cin >> x[i] >> y[i];
    
    std::vector<delaunay::Edge> all;
    for (size_t i=0; i<n; ++i) for (size_t j=i+1; j<n; ++j) all.push_back(delaunay::make_edge(i, j));
    
    delaunay::DelaunayTriangulation DT(x, y);
    DT.execute();
    std::vector<delaunay::Edge> edge = DT.get_edges();
    
    double mst_all = euclid_mst(x, y, all);
    double mst_DT = euclid_mst(x, y, edge);
    std::cout << "Using all edges: " << mst_all << std::endl;
    std::cout << "Using delaunay triangulation: " << mst_DT << std::endl;
    std::cout << "Error: " << std::abs(mst_all - mst_DT) << std::endl;
}
