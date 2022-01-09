#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

uint64_t mod_mul64(uint64_t a, uint64_t b, uint64_t mod) {
    assert(a < mod && b < mod);

    if (mod <= 1LLU << 32)
        return a * b % mod;

    if (mod <= 1LLU << 63) {
        uint64_t q = uint64_t((long double) a * b / mod);
        uint64_t result = a * b - q * mod;

        if (result > 1LLU << 63)
            result += mod;
        else if (result >= mod)
            result -= mod;

        return result;
    }

#ifdef __SIZEOF_INT128__
    return uint64_t(__uint128_t(a) * b % mod);
#endif

    assert(false);
}

uint64_t mod_pow64(uint64_t a, uint64_t b, uint64_t mod) {
    uint64_t result = 1;

    while (b > 0) {
        if (b & 1)
            result = mod_mul64(result, a, mod);

        a = mod_mul64(a, a, mod);
        b >>= 1;
    }

    return result;
}

bool miller_rabin(uint64_t n) {
    if (n < 2)
        return false;

    // Check small primes.
    for (uint64_t p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29})
        if (n % p == 0)
            return n == p;

    // https://miller-rabin.appspot.com/
    auto get_miller_rabin_bases = [&]() -> vector<uint64_t> {
        if (n < 341531) return {9345883071009581737LLU};
        if (n < 1050535501) return {336781006125, 9639812373923155};
        if (n < 350269456337) return {4230279247111683200, 14694767155120705706LLU, 16641139526367750375LLU};
        if (n < 55245642489451) return {2, 141889084524735, 1199124725622454117, 11096072698276303650LLU};
        if (n < 7999252175582851) return {2, 4130806001517, 149795463772692060, 186635894390467037, 3967304179347715805};
        if (n < 585226005592931977) return {2, 123635709730000, 9233062284813009, 43835965440333360, 761179012939631437, 1263739024124850375};
        return {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    };

    int r = __builtin_ctzll(n - 1);
    uint64_t d = (n - 1) >> r;

    for (uint64_t a : get_miller_rabin_bases()) {
        if (a % n == 0)
            continue;

        uint64_t x = mod_pow64(a % n, d, n);

        if (x == 1 || x == n - 1)
            continue;

        for (int i = 0; i < r - 1 && x != n - 1; i++)
            x = mod_mul64(x, x, n);

        if (x != n - 1)
            return false;
    }

    return true;
}

// Solution to https://www.spoj.com/problems/PON/
int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    uint64_t n;
    cin >> n;

    while (cin >> n)
        // cout << (miller_rabin(n) ? "YES" : "NO") << '\n';
        cout << n << ": " << (miller_rabin(n) ? "YES" : "NO") << '\n';
}
