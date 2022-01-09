// Solution to https://atcoder.jp/contests/dp/tasks/dp_v
// We need to paint every vertex of a tree white or black such that the black vertices are a single connected component.
// For every vertex v, find the number of ways to paint the tree such that v is painted black.
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

int N, MOD;
vector<vector<int>> adj;

// down[node] = number of ways to paint node's subtree where node must be painted black.
// up[node] = number of ways to color parent's subtree of node, including all white.
// combined[node] = answer for the whole tree when rooted at node.
vector<int64_t> down;
vector<int64_t> up;
vector<int64_t> combined;

void solve_down(int node, int parent) {
    vector<int> &children = adj[node];
    // Erase the edge from node to parent.
    children.erase(remove(children.begin(), children.end(), parent), children.end());
    down[node] = 1;

    for (int child : children) {
        assert(child != parent);
        solve_down(child, node);
        down[node] = down[node] * (down[child] + 1) % MOD;
    }
}

void solve_up(int node, int parent) {
    vector<int> &children = adj[node];
    int n = int(children.size());
    vector<int64_t> prefix(n + 1, 1);
    vector<int64_t> suffix(n + 1, 1);

    for (int i = 0; i < n; i++)
        prefix[i + 1] = prefix[i] * (down[children[i]] + 1) % MOD;

    for (int i = n - 1; i >= 0; i--)
        suffix[i] = suffix[i + 1] * (down[children[i]] + 1) % MOD;

    for (int i = 0; i < n; i++)
        up[children[i]] = (1 + prefix[i] * suffix[i + 1] % MOD * up[node]) % MOD;

    for (int child : children) {
        assert(child != parent);
        solve_up(child, node);
    }

    combined[node] = up[node] * down[node] % MOD;
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    cin >> N >> MOD;
    adj.assign(N, {});

    for (int i = 0; i < N - 1; i++) {
        int a, b;
        cin >> a >> b;
        a--; b--;
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    down.resize(N);
    solve_down(0, -1);

    up.resize(N);
    up[0] = 1;
    combined.resize(N);
    solve_up(0, -1);

    for (int i = 0; i < N; i++)
        cout << combined[i] << '\n';
}
