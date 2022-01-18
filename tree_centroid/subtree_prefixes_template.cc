// Solves the following: given a weighted tree of N nodes and a number K, how many paths in the tree have weight <= K?
// Runs in O(N log^2 N) time.
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

template<typename T>
struct fenwick_tree {
    static int highest_bit(unsigned x) {
        return x == 0 ? -1 : 31 - __builtin_clz(x);
    }

    int tree_n = 0;
    T tree_sum = T();
    vector<T> tree;

    fenwick_tree(int n = -1) {
        if (n >= 0)
            init(n);
    }

    void init(int n) {
        tree_n = n;
        tree_sum = T();
        tree.assign(tree_n + 1, T());
    }

    // O(n) initialization of the Fenwick tree.
    template<typename T_array>
    void build(const T_array &initial) {
        assert(int(initial.size()) == tree_n);
        tree_sum = T();

        for (int i = 1; i <= tree_n; i++) {
            tree[i] = initial[i - 1];
            tree_sum += initial[i - 1];

            for (int k = (i & -i) >> 1; k > 0; k >>= 1)
                tree[i] += tree[i - k];
        }
    }

    // index is in [0, tree_n).
    void update(int index, const T &change) {
        assert(0 <= index && index < tree_n);
        tree_sum += change;

        for (int i = index + 1; i <= tree_n; i += i & -i)
            tree[i] += change;
    }

    // Returns the sum of the range [0, count).
    T query(int count) const {
        count = min(count, tree_n);
        T sum = T();

        for (int i = count; i > 0; i -= i & -i)
            sum += tree[i];

        return sum;
    }

    // Returns the sum of the range [start, tree_n).
    T query_suffix(int start) const {
        return tree_sum - query(start);
    }

    // Returns the sum of the range [a, b).
    T query(int a, int b) const {
        return query(b) - query(a);
    }

    // Returns the element at index a in O(1) amortized across every index. Equivalent to query(a, a + 1).
    T get(int a) const {
        assert(0 <= a && a < tree_n);
        int above = a + 1;
        T sum = tree[above];
        above -= above & -above;

        while (a != above) {
            sum -= tree[a];
            a -= a & -a;
        }

        return sum;
    }

    bool set(int index, T value) {
        assert(0 <= index && index < tree_n);
        T current = get(index);

        if (current == value)
            return false;

        update(index, value - current);
        return true;
    }

    // Returns the largest p in `[0, tree_n]` such that `query(p) <= sum`. Returns -1 if no such p exists (`sum < 0`).
    // Can be used as an ordered set/multiset on indices in `[0, tree_n)` by using the tree as a 0/1 or frequency array:
    // `set(index, 1)` is equivalent to insert(index). `update(index, +1)` is equivalent to multiset.insert(index).
    // `set(index, 0)` or `update(index, -1)` are equivalent to erase(index).
    // `get(index)` provides whether index is present or not, or the count of index (if multiset).
    // `query(index)` provides the count of elements < index.
    // `find_last_prefix(k)` finds the k-th smallest element (0-indexed). Returns `tree_n` for `sum >= set.size()`.
    int find_last_prefix(T sum) const {
        if (sum < T())
            return -1;

        int prefix = 0;

        for (int k = highest_bit(tree_n); k >= 0; k--)
            if (prefix + (1 << k) <= tree_n && tree[prefix + (1 << k)] <= sum) {
                prefix += 1 << k;
                sum -= tree[prefix];
            }

        return prefix;
    }
};


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
    vector<int64_t> weighted_depth;
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
        weighted_depth.resize(N);
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
        weighted_depth[node] = weight;
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

    int64_t solve_subtree_prefixes() {
        int n = int(nodes.size());
        vector<int64_t> sorted_weights(n);

        for (int i = 0; i < n; i++)
            sorted_weights[i] = weighted_depth[nodes[i]];

        sort(sorted_weights.begin(), sorted_weights.end());
        fenwick_tree<int64_t> tree(n);
        int64_t pairs = 0;

        // Note that after DFS from the centroid, `nodes` is grouped by subtree.
        for (int i = 0, j = 0; i < n; i = j) {
            while (j < n && subroot[nodes[i]] == subroot[nodes[j]])
                j++;

            // Avoid paths within the same subtree by querying each subtree before updating it.
            for (int k = i; k < j; k++) {
                int64_t maximum = K - weighted_depth[nodes[k]];
                int position = int(upper_bound(sorted_weights.begin(), sorted_weights.end(), maximum) - sorted_weights.begin());
                pairs += tree.query(position);
            }

            for (int k = i; k < j; k++) {
                int position = int(lower_bound(sorted_weights.begin(), sorted_weights.end(), weighted_depth[nodes[k]]) - sorted_weights.begin());
                tree.update(position, +1);
            }
        }

        return pairs;
    }

    int64_t solve(int root) {
        root = centroid(root);

        for (int node : nodes)
            if (node != root)
                centroid_parent[node] = root;

        dfs(root);
        int64_t pairs = solve_subtree_prefixes();

        // TODO: determine whether the following code is needed to handle suffixes as well.

        // reverse(nodes.begin(), nodes.end());
        // pairs += solve_subtree_prefixes();

        // Recurse after solving root, so that edge erasures don't cause incorrect results.
        for (edge &e : adj[root]) {
            erase_edge(e.node, root);
            pairs += solve(e.node);
        }

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
