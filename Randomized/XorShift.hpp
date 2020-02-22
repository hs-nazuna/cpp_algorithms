#ifndef XOR_SHIFT_HPP
#define XOR_SHIFT_HPP

#include <cassert>

struct XorShift {
    unsigned int x, y, z, w;
    
    XorShift() { init(); }
    XorShift(unsigned int seed) { init(seed); }
    
    void init(unsigned int seed = 88675123) { x = 123456789; y = 362436069; z = 521288629; w = seed; }
    
    void setNext() {
        unsigned int t = x ^ (x << 11);
        x = y; y = z; z = w;
        w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }
    
    unsigned int nextInt() { setNext(); return w; }
    unsigned int nextInt(unsigned int n) { return nextInt() % n; }
    unsigned int nextInt(unsigned int a, unsigned int b) { assert(a <= b); return a + nextInt(b - a + 1); }
    double nextDouble() { return nextInt() * (1.0 / 0xFFFFFFFFu); }
};

#endif
