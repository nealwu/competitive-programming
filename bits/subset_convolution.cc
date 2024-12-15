// See https://codeforces.com/blog/entry/72488
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <vector>
using namespace std;

template<typename A, typename B> ostream& operator<<(ostream &os, const pair<A, B> &p) { return os << '(' << p.first << ", " << p.second << ')'; }
template<typename... Args> ostream& operator<<(ostream& os, const tuple<Args...>& t) { os << '('; apply([&os](const Args&... args) { size_t n = 0; ((os << args << (++n != sizeof...(Args) ? ", " : "")), ...); }, t); return os << ')'; }
template<typename T_container, typename T = typename enable_if<!is_same<T_container, string>::value, typename T_container::value_type>::type> ostream& operator<<(ostream &os, const T_container &v) { os << '{'; string sep; for (const T &x : v) os << sep << x, sep = ", "; return os << '}'; }

void dbg_out() { cerr << endl; }
template<typename Head, typename... Tail> void dbg_out(Head H, Tail... T) { cerr << ' ' << H; dbg_out(T...); }

#ifdef NEAL_DEBUG
#define dbg(...) cerr << '[' << __FILE__ << ':' << __LINE__ << "] (" << #__VA_ARGS__ << "):", dbg_out(__VA_ARGS__)
#else
#define dbg(...)
#endif

// For every mask, computes the sum of `values[sub]` where `sub` is a submask of mask.
template<typename T_out, typename T_in>
vector<T_out> submask_sums(const vector<T_in> &values) {
    int n = __builtin_ctz(int(values.size()));
    assert(int(values.size()) == 1 << n);
    vector<T_out> dp(values.begin(), values.end());

    // Broken profile DP where the intermediate DP state consists of the i-th suffix of the previous row and the i-th
    // prefix of the current row.
    for (int i = 0; i < n; i++)
        for (int base = 0; base < 1 << n; base += 1 << (i + 1))
            for (int mask = base; mask < base + (1 << i); mask++)
                dp[mask + (1 << i)] += dp[mask];

    return dp;
}

// Does the inverse of `submask_sums`; returns the input that produces the given output.
template<typename T_out, typename T_in>
vector<T_out> mobius_transform(const vector<T_in> &values) {
    int n = __builtin_ctz(int(values.size()));
    assert(int(values.size()) == 1 << n);
    vector<T_out> dp(values.begin(), values.end());

    for (int i = 0; i < n; i++)
        for (int base = 0; base < 1 << n; base += 1 << (i + 1))
            for (int mask = base; mask < base + (1 << i); mask++)
                dp[mask + (1 << i)] -= dp[mask];

    return dp;
}

template<typename F>
void iterate_bitmasks_with_popcount(int n, int k, F &&f) {
    if (k == 0) {
        f(0);
        return;
    }

    int mask = (1 << k) - 1;

    while (mask < 1 << n) {
        f(mask);
        int zeros = __builtin_ctz(mask);
        int ones = __builtin_ctz(~mask >> zeros);
        mask += (1 << zeros) + (1 << (ones - 1)) - 1;
    }
}

// Performs subset convolution, C[x | y] += A[x] * B[y] for all (x & y) = 0, in n^2 * 2^n time.
template<typename T_out, typename T_in>
vector<T_out> subset_convolution(const vector<T_in> &A, const vector<T_in> &B) {
    int n = __builtin_ctz(int(A.size()));
    assert(int(A.size()) == 1 << n && int(B.size()) == 1 << n);
    vector<vector<T_out>> FA(n + 1, vector<T_out>(1 << n, 0));
    vector<vector<T_out>> FB(n + 1, vector<T_out>(1 << n, 0));

    for (int mask = 0; mask < 1 << n; mask++) {
        FA[__builtin_popcount(mask)][mask] = A[mask];
        FB[__builtin_popcount(mask)][mask] = B[mask];
    }

    for (int c = 0; c < n; c++) {
        FA[c] = submask_sums<T_out>(FA[c]);
        FB[c] = submask_sums<T_out>(FB[c]);
    }

    vector<T_out> C(1 << n, 0);

    for (int c = 0; c <= n; c++) {
        vector<T_out> FC(1 << n, 0);

        // Add up all the ways to combine to c bits, with overlap.
        for (int i = 0; i <= c; i++)
            for (int mask = 0; mask < 1 << n; mask++)
                FC[mask] += FA[i][mask] * FB[c - i][mask];

        // Subtract out combinations that actually have fewer than c bits.
        if (c > 1)
            FC = mobius_transform<T_out>(FC);

        iterate_bitmasks_with_popcount(n, c, [&](int mask) -> void {
            C[mask] = FC[mask];
        });
    }

    return C;
}

// Performs reverse subset convolution, C[x] += A[x | y] * B[y] for all (x & y) = 0, in n^2 * 2^n time.
template<typename T_out, typename T_in>
vector<T_out> reverse_subset_convolution(vector<T_in> A, const vector<T_in> &B) {
    reverse(A.begin(), A.end());
    vector<T_out> result = subset_convolution<T_out>(A, B);
    reverse(result.begin(), result.end());
    return result;
}


template<bool add_one = false, typename T_vector>
void output_vector(const T_vector &v, int start = 0, int end = -1) {
    if (end < 0) end = int(v.size());

    for (int i = start; i < end; i++)
        if constexpr (add_one)
            cout << v[i] + 1 << (i < end - 1 ? ' ' : '\n');
        else
            cout << v[i] << (i < end - 1 ? ' ' : '\n');
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N;
    cin >> N;
    vector<int> A(1 << N), B(1 << N);

    for (auto &a : A)
        cin >> a;

    for (auto &b : B)
        cin >> b;

    output_vector(subset_convolution<int64_t>(A, B));
    output_vector(reverse_subset_convolution<int64_t>(A, B));
}
