#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

// Warning: this uses 195 MB for N = 500,000. If low on memory use the RMQ from memory_rmq_lca.cc or block_rmq_lca.cc
template<typename T, bool maximum_mode = false>
struct RMQ {
    static int highest_bit(unsigned x) {
        return x == 0 ? -1 : 31 - __builtin_clz(x);
    }

    int n = 0;
    vector<T> values;
    vector<vector<int>> range_low;

    RMQ(const vector<T> &_values = {}) {
        if (!_values.empty())
            build(_values);
    }

    // Note: when `values[a] == values[b]`, returns b.
    int better_index(int a, int b) const {
        return (maximum_mode ? values[b] < values[a] : values[a] < values[b]) ? a : b;
    }

    void build(const vector<T> &_values) {
        values = _values;
        n = int(values.size());
        int levels = highest_bit(n) + 1;
        range_low.resize(levels);

        for (int k = 0; k < levels; k++)
            range_low[k].resize(n - (1 << k) + 1);

        for (int i = 0; i < n; i++)
            range_low[0][i] = i;

        for (int k = 1; k < levels; k++)
            for (int i = 0; i <= n - (1 << k); i++)
                range_low[k][i] = better_index(range_low[k - 1][i], range_low[k - 1][i + (1 << (k - 1))]);
    }

    // Note: breaks ties by choosing the largest index.
    int query_index(int a, int b) const {
        assert(0 <= a && a < b && b <= n);
        int level = highest_bit(b - a);
        return better_index(range_low[level][a], range_low[level][b - (1 << level)]);
    }

    T query_value(int a, int b) const {
        return values[query_index(a, b)];
    }
};

template<typename T_weight>
struct weighted_LCA {
    struct edge {
        int node = -1;
        T_weight weight = 0;

        edge() {}

        edge(int _node, T_weight _weight) : node(_node), weight(_weight) {}
    };

    int n = 0;
    vector<vector<edge>> adj;
    vector<int> parent, depth, subtree_size;
    vector<T_weight> weighted_depth, up_weight;
    vector<int> euler, first_occurrence;
    vector<int> tour_start, tour_end;
    vector<int> tour_list, rev_tour_list;
    vector<int> heavy_root;
    vector<int> heavy_root_depth, heavy_root_parent;  // These two vectors serve purely to optimize get_kth_ancestor
    RMQ<int> rmq;
    bool built;

    weighted_LCA(int _n = 0) {
        init(_n);
    }

    // Warning: this does not call build().
    weighted_LCA(const vector<vector<edge>> &_adj) {
        init(_adj);
    }

    void init(int _n) {
        n = _n;
        adj.assign(n, {});
        parent.resize(n);
        depth.resize(n);
        subtree_size.resize(n);
        weighted_depth.resize(n);
        up_weight.assign(n, 0);
        first_occurrence.resize(n);
        tour_start.resize(n);
        tour_end.resize(n);
        tour_list.resize(n);
        heavy_root.resize(n);
        built = false;
    }

    // Warning: this does not call build().
    void init(const vector<vector<edge>> &_adj) {
        init(int(_adj.size()));
        adj = _adj;
    }

    void add_edge(int a, int b, T_weight weight) {
        adj[a].emplace_back(b, weight);
        adj[b].emplace_back(a, weight);
    }

    int degree(int v) const {
        return int(adj[v].size()) + (built && parent[v] >= 0);
    }

    void erase_edge(int from, int to) {
        for (edge &e : adj[from])
            if (e.node == to) {
                swap(e, adj[from].back());
                adj[from].pop_back();
                return;
            }
    }

    void dfs(int node, int par, T_weight weight) {
        parent[node] = par;
        depth[node] = par < 0 ? 0 : depth[par] + 1;
        subtree_size[node] = 1;
        weighted_depth[node] = weight;

        // Erase the edge to parent.
        erase_edge(node, par);

        for (edge &e : adj[node]) {
            up_weight[e.node] = e.weight;
            dfs(e.node, node, weight + e.weight);
            subtree_size[node] += subtree_size[e.node];
        }

        // Heavy-light subtree reordering.
        sort(adj[node].begin(), adj[node].end(), [&](const edge &a, const edge &b) {
            return subtree_size[a.node] > subtree_size[b.node];
        });
    }

    int tour;

    void tour_dfs(int node, bool heavy) {
        heavy_root[node] = heavy ? heavy_root[parent[node]] : node;
        first_occurrence[node] = int(euler.size());
        euler.push_back(node);
        tour_list[tour] = node;
        tour_start[node] = tour++;
        bool heavy_child = true;

        for (edge &e : adj[node]) {
            tour_dfs(e.node, heavy_child);
            euler.push_back(node);
            heavy_child = false;
        }

        tour_end[node] = tour;
    }

    void build(int root = -1, bool build_rmq = true) {
        parent.assign(n, -1);

        if (0 <= root && root < n)
            dfs(root, -1, 0);

        for (int i = 0; i < n; i++)
            if (i != root && parent[i] < 0)
                dfs(i, -1, 0);

        tour = 0;
        euler.clear();
        euler.reserve(2 * n);

        for (int i = 0; i < n; i++)
            if (parent[i] < 0) {
                tour_dfs(i, false);
                // Add a -1 in between connected components to help us detect when nodes aren't connected.
                euler.push_back(-1);
            }

        assert(int(euler.size()) == 2 * n);
        vector<int> euler_depth;
        euler_depth.reserve(euler.size());

        for (int node : euler)
            euler_depth.push_back(node < 0 ? node : depth[node]);

        if (build_rmq)
            rmq.build(euler_depth);

        euler_depth.clear();
        heavy_root_depth.resize(n);
        heavy_root_parent.resize(n);

        for (int i = 0; i < n; i++) {
            heavy_root_depth[i] = depth[heavy_root[i]];
            heavy_root_parent[i] = parent[heavy_root[i]];
        }

        rev_tour_list = tour_list;
        reverse(rev_tour_list.begin(), rev_tour_list.end());
        built = true;
    }

    pair<int, array<int, 2>> get_diameter() const {
        assert(built);

        // We find the maximum of depth[u] - 2 * depth[x] + depth[v] where u, x, v occur in order in the Euler tour.
        pair<int, int> u_max = {-1, -1};
        pair<int, int> ux_max = {-1, -1};
        pair<int, array<int, 2>> uxv_max = {-1, {-1, -1}};

        for (int node : euler) {
            if (node < 0) break;
            u_max = max(u_max, {depth[node], node});
            ux_max = max(ux_max, {u_max.first - 2 * depth[node], u_max.second});
            uxv_max = max(uxv_max, {ux_max.first + depth[node], {ux_max.second, node}});
        }

        return uxv_max;
    }

    // Returns the center(s) of the tree (the midpoint(s) of the diameter).
    array<int, 2> get_center() const {
        pair<int, array<int, 2>> diam = get_diameter();
        int length = diam.first, a = diam.second[0], b = diam.second[1];
        return {get_kth_node_on_path(a, b, length / 2), get_kth_node_on_path(a, b, (length + 1) / 2)};
    }

    pair<T_weight, array<int, 2>> get_weighted_diameter() const {
        assert(built);

        // We find the maximum of weighted_depth[u] - 2 * weighted_depth[x] + weighted_depth[v]
        // where u, x, v occur in order in the Euler tour.
        pair<T_weight, int> u_max = {-1, -1};
        pair<T_weight, int> ux_max = {-1, -1};
        pair<T_weight, array<int, 2>> uxv_max = {-1, {-1, -1}};

        for (int node : euler) {
            if (node < 0) break;
            u_max = max(u_max, {weighted_depth[node], node});
            ux_max = max(ux_max, {u_max.first - 2 * weighted_depth[node], u_max.second});
            uxv_max = max(uxv_max, {ux_max.first + weighted_depth[node], {ux_max.second, node}});
        }

        return uxv_max;
    }

    // Note: returns -1 if `a` and `b` aren't connected.
    int get_lca(int a, int b) const {
        a = first_occurrence[a];
        b = first_occurrence[b];

        if (a > b)
            swap(a, b);

        return euler[rmq.query_index(a, b + 1)];
    }

    bool is_ancestor(int a, int b) const {
        return tour_start[a] <= tour_start[b] && tour_start[b] < tour_end[a];
    }

    bool on_path(int x, int a, int b) const {
        return (is_ancestor(x, a) || is_ancestor(x, b)) && is_ancestor(get_lca(a, b), x);
    }

    int get_dist(int a, int b) const {
        return depth[a] + depth[b] - 2 * depth[get_lca(a, b)];
    }

    T_weight get_weighted_dist(int a, int b) const {
        return weighted_depth[a] + weighted_depth[b] - 2 * weighted_depth[get_lca(a, b)];
    }

    // Returns the child of `a` that is an ancestor of `b`. Assumes `a` is a strict ancestor of `b`.
    int child_ancestor(int a, int b) const {
        assert(a != b && is_ancestor(a, b));

        // Note: this depends on RMQ breaking ties by latest index.
        int child = euler[rmq.query_index(first_occurrence[a], first_occurrence[b] + 1) + 1];
        // assert(parent[child] == a && is_ancestor(child, b));
        return child;
    }

    int get_kth_ancestor(int a, int k) const {
        if (k > depth[a])
            return -1;

        int goal = depth[a] - k;

        while (heavy_root_depth[a] > goal)
            a = heavy_root_parent[a];

        return tour_list[tour_start[a] + goal - depth[a]];
    }

    int get_kth_node_on_path(int a, int b, int k) const {
        int anc = get_lca(a, b);
        int first_half = depth[a] - depth[anc];
        int second_half = depth[b] - depth[anc];

        if (k < 0 || k > first_half + second_half)
            return -1;

        if (k < first_half)
            return get_kth_ancestor(a, k);
        else
            return get_kth_ancestor(b, first_half + second_half - k);
    }

    // Note: this is the LCA of any two nodes out of three when the third node is the root.
    // It is also the node with the minimum sum of distances to all three nodes (the centroid of the three nodes).
    int get_common_node(int a, int b, int c) const {
        // Return the deepest node among lca(a, b), lca(b, c), and lca(c, a).
        int x = get_lca(a, b);
        int y = get_lca(b, c);
        int z = get_lca(c, a);
        return x ^ y ^ z;
    }

    // Given a subset of k tree nodes, computes the minimal subtree that contains all the nodes (at most 2k - 1 nodes).
    // Returns a list of {node, parent} for every node in the subtree sorted by tour index. Runs in O(k log k).
    // Note that all parents also appear as a node in the return value, and nodes[0].first is the compressed root.
    vector<pair<int, int>> compress_tree(vector<int> nodes) const {
        if (nodes.empty())
            return {};

        auto &&compare_tour = [&](int a, int b) { return tour_start[a] < tour_start[b]; };
        sort(nodes.begin(), nodes.end(), compare_tour);
        int k = int(nodes.size());

        for (int i = 0; i < k - 1; i++)
            nodes.push_back(get_lca(nodes[i], nodes[i + 1]));

        sort(nodes.begin() + k, nodes.end(), compare_tour);
        inplace_merge(nodes.begin(), nodes.begin() + k, nodes.end(), compare_tour);
        nodes.erase(unique(nodes.begin(), nodes.end()), nodes.end());
        vector<pair<int, int>> result = {{nodes[0], -1}};

        for (int i = 1; i < int(nodes.size()); i++)
            result.emplace_back(nodes[i], get_lca(nodes[i], nodes[i - 1]));

        return result;
    }
};


int main() {
    int N, Q;
    scanf("%d %d", &N, &Q);
    weighted_LCA<int64_t> lca(N);

    for (int i = 0; i < N - 1; i++) {
        int a, b;
        long long weight;
        scanf("%d %d %lld", &a, &b, &weight);
        a--; b--;
        lca.add_edge(a, b, weight);
    }

    lca.build();
    auto diameter = lca.get_weighted_diameter();
    printf("%lld\n", (long long) diameter.first);
    assert(diameter.first == lca.get_weighted_dist(diameter.second[0], diameter.second[1]));

    for (int q = 0; q < Q; q++) {
        int a, b;
        scanf("%d %d", &a, &b);
        a--; b--;
        printf("%lld\n", (long long) lca.get_weighted_dist(a, b));
    }
}
