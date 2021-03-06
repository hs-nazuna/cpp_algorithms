#ifndef TSUKASA_DIARY_COMBINATION_HPP
#define TSUKASA_DIARY_COMBINATION_HPP

#include <vector>

using lint = long long int;

class Combination {
private:
    const int N;
    const lint MOD;
    std::vector<lint> fact, inv;
    
    lint powmod(lint x, lint k) {
        lint res = 1;
        while (k) {
            if (k & 1) res = (res * x) % MOD;
            x = (x * x) % MOD;
            k >>= 1;
        }
        return res;
    }
    
    void precomputation() {
        fact[0] = fact[1] = 1;
        for (int i=2; i<=N; ++i) fact[i] = (fact[i - 1] * i) % MOD;
        for (int i=0; i<=N; ++i) inv[i] = powmod(fact[i], MOD - 2);
    }
    
public:
    Combination(int N, lint MOD) : N(N), MOD(MOD), fact(N + 1), inv(N + 1) {
        precomputation();
    }

    lint P(int n, int k) {
        if (k == 0) return 1;
        if (n < k) return 0;
        return (fact[n] * inv[n - k]) % MOD;
    }
    
    lint C(int n, int k) {
        if (k == 0) return 1;
        if (n < k) return 0;
        lint res = fact[n];
        res = (res * inv[n - k]) % MOD;
        res = (res * inv[k]) % MOD;
        return res;
    }

    lint H(int n, int k) {
        return (n + k - 1, k);
    }
};

#endif
