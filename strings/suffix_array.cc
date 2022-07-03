#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;

template<typename T, bool maximum_mode = false>
struct RMQ {
    static int highest_bit(unsigned x) {
        return x == 0 ? -1 : 31 - __builtin_clz(x);
    }

    int n = 0;
    vector<vector<T>> range_min;

    RMQ(const vector<T> &values = {}) {
        if (!values.empty())
            build(values);
    }

    static T better(T a, T b) {
        return maximum_mode ? max(a, b) : min(a, b);
    }

    void build(const vector<T> &values) {
        n = int(values.size());
        int levels = highest_bit(n) + 1;
        range_min.resize(levels);

        for (int k = 0; k < levels; k++)
            range_min[k].resize(n - (1 << k) + 1);

        if (n > 0)
            range_min[0] = values;

        for (int k = 1; k < levels; k++)
            for (int i = 0; i <= n - (1 << k); i++)
                range_min[k][i] = better(range_min[k - 1][i], range_min[k - 1][i + (1 << (k - 1))]);
    }

    T query_value(int a, int b) const {
        assert(0 <= a && a < b && b <= n);
        int level = highest_bit(b - a);
        return better(range_min[level][a], range_min[level][b - (1 << level)]);
    }
};

template<typename T_string = string>
struct suffix_array {
    int n = 0;
    T_string str;
    vector<int> suffix;
    vector<int> rank;
    vector<int> lcp;
    RMQ<int> rmq;

    suffix_array() {}

    suffix_array(const T_string &_str, bool build_rmq = true) {
        build(_str, build_rmq);
    }

    void build(const T_string &_str, bool build_rmq = true) {
        str = _str;
        n = int(str.size());
        suffix.resize(n);
        iota(suffix.begin(), suffix.end(), 0);
        bool large_alphabet = false;

        for (int i = 0; i < n; i++)
            if (str[i] < 0 || str[i] >= 128)
                large_alphabet = true;

        // Sort each suffix by the first character.
        if (large_alphabet) {
            sort(suffix.begin(), suffix.end(), [&](int a, int b) {
                return str[a] < str[b];
            });
        } else {
            vector<int> freq(128, 0);

            for (int i = 0; i < n; i++)
                freq[str[i]]++;

            for (int c = 1; c < 128; c++)
                freq[c] += freq[c - 1];

            for (int i = 0; i < n; i++)
                suffix[--freq[str[i]]] = i;
        }

        // Compute the rank of each suffix. Tied suffixes share the same rank.
        rank.resize(n);
        rank[suffix[0]] = 0;

        for (int i = 1; i < n; i++)
            rank[suffix[i]] = str[suffix[i]] == str[suffix[i - 1]] ? rank[suffix[i - 1]] : i;

        bool done = false;

        for (int len = 1; len < n && !done; len *= 2) {
            // next_index[i] = the next index to use for a suffix of rank i. We insert them in order of the rank of the
            // suffix that comes len characters after the current suffix.
            vector<int> next_index(n);
            iota(next_index.begin(), next_index.end(), 0);

            // Compute the suffix array for 2 * len. Suffixes of length <= len are prioritized first.
            vector<int> nsuffix(n);

            for (int i = n - len; i < n; i++)
                nsuffix[next_index[rank[i]]++] = i;

            for (int i = 0; i < n; i++) {
                int before = suffix[i] - len;

                if (before >= 0)
                    nsuffix[next_index[rank[before]]++] = before;
            }

            suffix.swap(nsuffix);
            nsuffix.clear();

            // Compute the rank array for 2 * len.
            vector<int> nrank(n);
            nrank[suffix[0]] = 0;
            done = true;

            for (int i = 1; i < n; i++) {
                int s = suffix[i], before = suffix[i - 1];

                if (max(s, before) + len < n && rank[s] == rank[before] && rank[s + len] == rank[before + len]) {
                    nrank[s] = nrank[before];
                    done = false;
                } else {
                    nrank[s] = i;
                }
            }

            rank.swap(nrank);
        }

        compute_lcp();

        if (build_rmq)
            rmq.build(lcp);
    }

    void compute_lcp() {
        lcp.assign(n, 0);
        int match = 0;

        for (int i = 0; i < n; i++) {
            if (rank[i] == 0)
                continue;

            int a = suffix[rank[i]] + match;
            int b = suffix[rank[i] - 1] + match;

            while (a < n && b < n && str[a++] == str[b++])
                match++;

            // lcp[r] = the longest common prefix length of the suffixes starting at suffix[r] and at suffix[r - 1].
            // Note that lcp[0] is always 0.
            lcp[rank[i]] = match;
            match = max(match - 1, 0);
        }
    }

    int get_lcp_from_ranks(int a, int b) const {
        if (a == b)
            return n - suffix[a];

        if (a > b)
            swap(a, b);

        return rmq.query_value(a + 1, b + 1);
    }

    int get_lcp(int a, int b) const {
        if (a >= n || b >= n)
            return 0;

        if (a == b)
            return n - a;

        return get_lcp_from_ranks(rank[a], rank[b]);
    }

    // Compares the substrings starting at `a` and `b` up to `length` in O(1).
    int compare(int a, int b, int length = -1) const {
        if (length < 0)
            length = n;

        if (a == b)
            return 0;

        int common = get_lcp(a, b);

        if (common >= length)
            return 0;

        if (a + common >= n || b + common >= n)
            return a + common >= n ? -1 : 1;

        return str[a + common] < str[b + common] ? -1 : (str[a + common] == str[b + common] ? 0 : 1);
    }
};


#include <numeric>

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    string task;
    cin >> task;

    string str;
    cin >> str;
    int N = int(str.size());
    suffix_array<string> SA(str);

    if (task == "suffix_array") {
        for (int i = 0; i < N; i++)
            cout << SA.suffix[i] << (i < N - 1 ? ' ' : '\n');

        for (int i = 0; i < N; i++)
            cout << SA.lcp[i] << (i < N - 1 ? ' ' : '\n');
    } else if (task == "distinct_substrings") {
        cout << int64_t(N) * (N + 1) / 2 - accumulate(SA.lcp.begin(), SA.lcp.end(), 0LL) << '\n';
    } else {
        assert(false);
    }
}
