#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
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

// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2016/p0200r0.html
template<class Fun> class y_combinator_result {
    Fun fun_;
public:
    template<class T> explicit y_combinator_result(T &&fun): fun_(std::forward<T>(fun)) {}
    template<class ...Args> decltype(auto) operator()(Args &&...args) { return fun_(std::ref(*this), std::forward<Args>(args)...); }
};
template<class Fun> decltype(auto) y_combinator(Fun &&fun) { return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun)); }


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

// Compresses the values in arr to be in the range [0, n).
template<typename T_out = int, typename T>
pair<vector<T_out>, vector<T>> compress_array(const vector<T> &arr) {
    int n = int(arr.size());
    vector<pair<T, int>> sorted(n);

    for (int i = 0; i < n; i++)
        sorted[i] = {arr[i], i};

    sort(sorted.begin(), sorted.end(), [](const pair<T, int> &x, const pair<T, int> &y) -> bool {
        return x.first < y.first;
    });

    vector<T_out> compressed(n);
    vector<T> sorted_values;
    sorted_values.reserve(n);
    int current = 0;

    for (int i = 0, j = 0; i < n; i = j) {
        while (j < n && sorted[i].first == sorted[j].first)
            compressed[sorted[j++].second] = current;

        sorted_values.push_back(sorted[i].first);
        current++;
    }

    return {compressed, sorted_values};
}


// Compresses the values in arr to be in the range [0, n).
template<typename T_out = int, typename T>
pair<vector<T_out>, vector<T>> compress_array_binary_search(const vector<T> &arr) {
    int n = int(arr.size());
    vector<T> sorted = arr;
    sort(sorted.begin(), sorted.end());
    sorted.erase(unique(sorted.begin(), sorted.end()), sorted.end());
    vector<T_out> compressed(n);

    for (int i = 0; i < n; i++)
        compressed[i] = int(lower_bound(sorted.begin(), sorted.end(), arr[i]) - sorted.begin());

    return {compressed, sorted};
}


uint64_t random_address() { char *p = new char; delete p; return uint64_t(p); }

const uint64_t SEED = chrono::steady_clock::now().time_since_epoch().count() * (random_address() | 1);
mt19937_64 rng(SEED);

template<typename T_out, typename T>
uint64_t compute_hash(const pair<vector<T_out>, vector<T>> &result) {
    uint64_t hash = 0;

    for (const T_out &x : result.first)
        hash = 123456789 * hash + x;

    for (const T &x : result.second)
        hash = 123456789 * hash + x;

    return hash;
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    cerr << fixed << setprecision(3);

    int N;
    int64_t A_MAX;

    if (argc > 2) {
        N = stoi(argv[1]);
        A_MAX = stoll(argv[2]);
    } else {
        N = int(1e6);
        A_MAX = int64_t(1e10);
    }

    vector<int64_t> A(N);

    for (auto &a : A)
        a = rng() % A_MAX;

    long double begin;
    pair<vector<unsigned>, vector<int64_t>> result;

    begin = clock();
    result = compress_array<unsigned>(A);
    cerr << (clock() - begin) / CLOCKS_PER_SEC << 's' << endl;
    uint64_t hash0 = compute_hash(result);
    cerr << "hash = " << hash0 << endl;

    begin = clock();
    result = compress_array_binary_search<unsigned>(A);
    cerr << (clock() - begin) / CLOCKS_PER_SEC << 's' << endl;
    uint64_t hash1 = compute_hash(result);
    cerr << "hash = " << hash1 << endl;

    assert(hash0 == hash1);
}
