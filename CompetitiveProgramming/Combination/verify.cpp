// verify @ https://atcoder.jp/contests/abc156/submissions/10293759

#include <bits/stdc++.h>
using namespace std;

using lint = long long int;

class Combination {
private:
    const int N;
    const lint MOD;
    std::vector< lint > fact, inv;
    
    lint powmod(lint x, lint k) {
        lint res = 1;
        while (k) {
            if (k&1) res = (res * x) % MOD;
            x = (x * x) % MOD;
            k >>= 1;
        }
        return res;
    }
    
    void precomputation() {
        fact[0] = fact[1] = 1;
        for (int i=2; i<=N; ++i) fact[i] = (fact[i-1] * i) % MOD;
        for (int i=0; i<=N; ++i) inv[i] = powmod(fact[i], MOD-2);
    }
    
public:
    Combination(int N_, lint MOD_) : N(N_), MOD(MOD_), fact(N_), inv(N_) {
        precomputation();
    }
    
    lint compute(int n, int k, bool duplicate = false) {
        if (duplicate) return compute(n+k-1, k);
        if (n < k) return 0;
        lint res = fact[n];
        res = (res * inv[n-k]) % MOD;
        res = (res * inv[k]) % MOD;
        return res;
    }
};

int main() {
    const lint MOD = (lint)1e9 + 7;

    int n, k;
    cin >> n >> k;

    Combination comb(2 * n, MOD);

    lint ans = 1;
    for (int i=1; i<=min(n, k); ++i) {
        lint add = comb.compute(n, i);
        add = (add * comb.compute(n-i, i, true)) % MOD;
        ans = (ans + add) % MOD;
    }
    cout << ans << endl;
}
