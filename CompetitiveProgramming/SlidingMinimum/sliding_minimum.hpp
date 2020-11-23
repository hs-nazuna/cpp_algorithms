#ifndef TSUKASA_DIARY_SLIDING_MINIMUM_HPP
#define TSUKASA_DIARY_SLIDING_MINIMUM_HPP

#include <vector>
#include <deque>

template<typename T>
std::vector<int> sliding_minimum(const std::vector<T>& tar, int K) {
    std::vector<int> ret;    
    std::deque<int> deq;
    for (int i=0, n=tar.size(); i<n; ++i) {
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
