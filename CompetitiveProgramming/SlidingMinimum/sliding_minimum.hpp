#ifndef TSUKASA_DIARY_SLIDING_MINIMUM_HPP
#define TSUKASA_DIARY_SLIDING_MINIMUM_HPP

#include <vector>
#include <deque>
#include <cassert>

template<typename T>
std::vector<size_t> sliding_minimum(const std::vector<T>& tar, size_t K) {
    assert(K > 0);
    
    std::vector<size_t> ret;
    
    std::deque<size_t> deq;
    for (size_t i=0, n=tar.size(); i<n; ++i) {
        while (!deq.empty() and tar[deq.back()] >= tar[i]) deq.pop_back();
        deq.push_back(i);
        if (i >= K - 1) {
            ret.push_back(deq.front());
            if (deq.front() + K - 1 == i) deq.pop_front();
        }
    }
    
    return ret;
}

#endif
