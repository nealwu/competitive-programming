// Solution to https://codeforces.com/contest/1249/problem/F
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

template<typename T1, typename T2>
bool maximize(T1 &a, const T2 &b) {
    if (a < b) {
        a = b;
        return true;
    }

    return false;
}

const int INF = int(1e9) + 5;

struct result {
    // dp[d] = the maximum weight to choose some nodes out of our subtree such that the highest node (closest to the
    // root) is depth d.
    vector<int> dp;

    result(int value = -INF) {
        dp = {value};
    }

    int size() const {
        return int(dp.size());
    }

    int maximum() const {
        return *max_element(dp.begin(), dp.end());
    }
};

int N, K;
vector<int> A;
vector<vector<int>> adj;
vector<result> results;

void extend(result &root) {
    // Shift the d dimension up by one.
    root.dp.insert(root.dp.begin(), -INF);
}

void attach(result &root, result &child) {
    result combined;
    combined.dp.assign(max(root.size(), child.size()), -INF);

    for (int i = 0; i < root.size(); i++)
        maximize(combined.dp[i], root.dp[i]);

    for (int i = 0; i < child.size(); i++)
        maximize(combined.dp[i], child.dp[i]);

    for (int x = 0; x < root.size(); x++)
        for (int y = max(K - x, 0); y < child.size(); y++)
            maximize(combined.dp[min(x, y)], root.dp[x] + child.dp[y]);

    root = combined;
}

void dfs(int node, int parent) {
    result &current = results[node];
    current = result(A[node]);

    for (int neigh : adj[node])
        if (neigh != parent) {
            dfs(neigh, node);
            result &child = results[neigh];
            extend(child);
            attach(current, child);
        }
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    cin >> N >> K;
    K++;
    A.resize(N);
    adj.assign(N, {});

    for (auto &a : A)
        cin >> a;

    for (int i = 0; i < N - 1; i++) {
        int u, v;
        cin >> u >> v;
        u--; v--;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    results.assign(N, {});
    dfs(0, -1);
    cout << results[0].maximum() << '\n';
}
