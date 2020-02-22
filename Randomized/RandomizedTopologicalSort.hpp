#ifndef RANDOMIZED_TOPOLOGICAL_SORT_HPP
#define RANDOMIZED_TOPOLOGICAL_SORT_HPP

#include <vector>
#include "UniformRandomSelector.hpp"

class RandomizedTopologicalSort {
private:
    URS< int > urs_que;
    
public:
    RandomizedTopologicalSort() {}
    
    std::vector< int > generate(std::vector< std::vector< int > >& adj_list) {
        urs_que.clean();
        
        int n = adj_list.size();
        std::vector< int > deg(n, 0);
        for (int v=0; v<n; ++v) for (int to : adj_list[v]) ++deg[to];
        for (int v=0; v<n; ++v) if (deg[v] == 0) urs_que.add(v);
        
        int i = 0;
        std::vector< int > ret(n, -1);
        while (urs_que.size() > 0) {
            int v = urs_que.select();
            ret[i++] = v;
            for (int to : adj_list[v]) {
                --deg[to];
                if (deg[to] == 0) urs_que.add(to);
            }
        }		
        return ret;
    }
};

#endif
