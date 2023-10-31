// Solution to https://atcoder.jp/contests/dp/tasks/dp_p
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

const int MOD = int(1e9) + 7;

struct result {
    int64_t black = 1, white = 1;

    int64_t either() const {
        return (black + white) % MOD;
    }
};

struct tree_dp {
    int N;
    vector<vector<int>> adj;

    tree_dp(int _N = 0) {
        init(_N);
    }

    void init(int _N) {
        N = _N;
        adj.assign(N, {});
    }

    void add_edge(int u, int v) {
        adj[u].emplace_back(v);
        adj[v].emplace_back(u);
    }

    // Warning: this setup is only recommended when the `result` struct is O(1) space.
    // Otherwise there is too much copying of `result` structs. In that case, use arrays_template.cc.
    result extend(result root) {
        return root;
    }

    result attach(result root, result child) {
        result combined;
        combined.black = root.black * child.white % MOD;
        combined.white = root.white * child.either() % MOD;
        return combined;
    }

    result dfs(int node, int parent) {
        result current;

        for (int neighbor : adj[node])
            if (neighbor != parent) {
                result child = dfs(neighbor, node);
                current = attach(current, extend(child));
            }

        return current;
    }

    result solve() {
        return dfs(0, -1);
    }
};

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N;
    cin >> N;
    tree_dp dp(N);

    for (int i = 0; i < N - 1; i++) {
        int u, v;
        cin >> u >> v;
        u--; v--;
        dp.add_edge(u, v);
    }

    cout << dp.solve().either() << '\n';
}
