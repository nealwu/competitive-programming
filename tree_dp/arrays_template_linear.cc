// Solution to https://codeforces.com/contest/1249/problem/F
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <queue>
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

struct result {
    // dp[d] = the maximum possible weight when choosing some nodes out of our subtree such that the closest node to the
    // subtree root has depth at least d. In particular, dp[d] >= dp[d + 1] for all d.
    deque<int> dp;

    result(int value = 0) {
        dp = {value};
    }

    int size() const {
        return int(dp.size());
    }

    int get(int index) const {
        return index < size() ? dp[index] : 0;
    }
};

int N, K;
vector<int> A;
vector<vector<int>> adj;
vector<result> results;

void extend(result &root) {
    // Shift the d dimension up by one.
    root.dp.push_front(root.dp.front());
}

void attach(result &root, result &child) {
    if (root.size() < child.size())
        swap(root, child);

    vector<int> combined_dp(child.dp.begin(), child.dp.end());

    for (int d = 0; d < child.size(); d++) {
        int min_other = max(K - d, d);
        maximize(combined_dp[d], root.dp[d] + child.get(min_other));
        maximize(combined_dp[d], root.get(min_other) + child.dp[d]);
    }

    int maximum = 0;

    for (int i = child.size() - 1; i >= 0; i--) {
        maximize(maximum, combined_dp[i]);
        maximize(root.dp[i], maximum);
    }
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
    cout << results[0].dp.front() << '\n';
}
