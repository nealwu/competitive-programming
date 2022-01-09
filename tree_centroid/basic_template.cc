#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

struct edge {
    int node;
    int64_t weight;

    edge(int _node = -1, int64_t _weight = 0) : node(_node), weight(_weight) {}
};

struct centroid_decomposition {
    int N;
    vector<vector<edge>> adj;
    vector<int> depth;
    vector<int> subtree_size;
    vector<int> centroid_parent;
    vector<int> subroot;
    vector<int> nodes;

    centroid_decomposition(int _N = 0) {
        init(_N);
    }

    void init(int _N) {
        N = _N;
        adj.assign(N, {});
        depth.resize(N);
        subtree_size.resize(N);
        centroid_parent.assign(N, -1);
        subroot.resize(N);
    }

    void add_edge(int u, int v, int64_t weight = 0) {
        assert(u != v);
        adj[u].emplace_back(v, weight);
        adj[v].emplace_back(u, weight);
    }

    void erase_edge(int from, int to) {
        for (edge &e : adj[from])
            if (e.node == to) {
                swap(e, adj[from].back());
                adj[from].pop_back();
                return;
            }

        assert(false);
    }

    int dfs(int node, int parent = -1, int sub = -1, int64_t weight = 0) {
        if (parent < 0) {
            sub = node;
            nodes.clear();
        }

        depth[node] = parent < 0 ? 0 : depth[parent] + 1;
        subtree_size[node] = 1;
        subroot[node] = sub;
        nodes.push_back(node);

        for (edge &e : adj[node])
            if (e.node != parent)
                subtree_size[node] += dfs(e.node, node, parent < 0 ? e.node : sub, weight + e.weight);

        return subtree_size[node];
    }

    int centroid(int root) {
        int n = dfs(root);
        bool found;

        // Repeatedly move to the subtree that is at least half of the tree, if such a subtree exists.
        do {
            found = false;

            for (edge &e : adj[root])
                if (subtree_size[e.node] < subtree_size[root] && 2 * subtree_size[e.node] >= n) {
                    root = e.node;
                    found = true;
                    break;
                }
        } while (found);

        return root;
    }

    void solve(int root) {
        root = centroid(root);

        for (int node : nodes)
            if (node != root)
                centroid_parent[node] = root;

        // TODO: either compute the answer for the whole tree here by calling `dfs(root)`, or compute it for each
        // subtree below by calling `dfs(e.node)`.

        for (edge &e : adj[root]) {
            erase_edge(e.node, root);
            // TODO: either compute the answer for the subtree of `e.node` here, or compute it for the whole tree above.
            // If computing for the whole tree above, we can move the `erase_edge` call to the loop below instead.
        }

        // Recurse after solving root, so that edge erasures don't cause incorrect results.
        for (edge &e : adj[root])
            solve(e.node);
    }
};

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N;
    cin >> N;
    centroid_decomposition CD(N);

    for (int i = 0; i < N - 1; i++) {
        int u, v;
        cin >> u >> v;
        u--; v--;
        CD.add_edge(u, v);
    }

    CD.solve(0);
}
