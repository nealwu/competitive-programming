// Solves the following: given a weighted tree of N nodes and a number K, how many paths in the tree have weight <= K?
// Runs in O(N log^2 N) time.
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
    int64_t K;
    vector<vector<edge>> adj;
    vector<int> depth;
    vector<int> subtree_size;
    vector<int> centroid_parent;
    vector<int> subroot;
    vector<int> nodes;
    vector<int64_t> weights;

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
            weights.clear();
        }

        depth[node] = parent < 0 ? 0 : depth[parent] + 1;
        subtree_size[node] = 1;
        subroot[node] = sub;
        nodes.push_back(node);
        weights.push_back(weight);

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

    int64_t count_pairs(int root, int64_t weight_max) {
        int n = dfs(root);
        int64_t pairs = 0;
        sort(weights.begin(), weights.end());
        assert(int(weights.size()) == n);

        for (int i = 0, j = n - 1; i < j; i++) {
            while (j > i && weights[i] + weights[j] > weight_max)
                j--;

            pairs += j - i;
        }

        return pairs;
    }

    int64_t solve(int root) {
        root = centroid(root);

        for (int node : nodes)
            if (node != root)
                centroid_parent[node] = root;

        // Compute the crossing pairs by counting all pairs and then subtracting pairs within the same subtree.
        int64_t pairs = count_pairs(root, K);

        for (edge &e : adj[root]) {
            erase_edge(e.node, root);
            pairs -= count_pairs(e.node, K - 2 * e.weight);
        }

        // Recurse after solving root, so that edge erasures don't cause incorrect results.
        for (edge &e : adj[root])
            pairs += solve(e.node);

        return pairs;
    }
};

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N;
    int64_t K;
    cin >> N >> K;
    centroid_decomposition CD(N);
    CD.K = K;

    for (int i = 0; i < N - 1; i++) {
        int u, v;
        int64_t weight;
        cin >> u >> v >> weight;
        u--; v--;
        CD.add_edge(u, v, weight);
    }

    cout << CD.solve(0) << '\n';
}
