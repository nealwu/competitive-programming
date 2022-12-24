#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

template<typename A, typename B> ostream& operator<<(ostream &os, const pair<A, B> &p) { return os << '(' << p.first << ", " << p.second << ')'; }
template<typename T_container, typename T = typename enable_if<!is_same<T_container, string>::value, typename T_container::value_type>::type> ostream& operator<<(ostream &os, const T_container &v) { os << '{'; string sep; for (const T &x : v) os << sep << x, sep = ", "; return os << '}'; }

void dbg_out() { cerr << endl; }
template<typename Head, typename... Tail> void dbg_out(Head H, Tail... T) { cerr << ' ' << H; dbg_out(T...); }

#ifdef NEAL_DEBUG
#define dbg(...) cerr << '[' << __FILE__ << ':' << __LINE__ << "] (" << #__VA_ARGS__ << "):", dbg_out(__VA_ARGS__)
#else
#define dbg(...)
#endif

struct biconnected_components {
    struct edge {
        int node, index;

        edge() {}

        edge(int _node, int _index) : node(_node), index(_index) {}
    };

    int n, edges;
    vector<vector<edge>> adj;
    vector<array<int, 2>> edge_list;
    vector<int> tour_start;
    vector<int> low_link;

    vector<bool> visited;
    vector<bool> is_cut;
    vector<bool> is_bridge;
    vector<int> stack;
    vector<vector<int>> components;
    int tour;

    biconnected_components(int _n = 0) {
        init(_n);
    }

    void init(int _n) {
        n = _n;
        edges = 0;
        adj.assign(n, {});
        edge_list.clear();
        tour_start.resize(n);
        low_link.resize(n);
    }

    void add_edge(int a, int b) {
        adj[a].emplace_back(b, edges);
        adj[b].emplace_back(a, edges);
        edge_list.push_back({a, b});
        edges++;
    }

    void dfs(int node, int parent) {
        assert(!visited[node]);
        visited[node] = true;
        tour_start[node] = tour++;
        low_link[node] = tour_start[node];
        is_cut[node] = false;
        int parent_count = 0, children = 0;

        for (edge &e : adj[node]) {
            // Skip the first edge to the parent, but allow multi-edges.
            if (e.node == parent && parent_count++ == 0)
                continue;

            if (visited[e.node]) {
                // e.node is a candidate for low_link.
                low_link[node] = min(low_link[node], tour_start[e.node]);

                if (tour_start[e.node] < tour_start[node])
                    stack.push_back(node);
            } else {
                int size = int(stack.size());
                dfs(e.node, node);
                children++;

                // e.node is part of our subtree.
                low_link[node] = min(low_link[node], low_link[e.node]);

                if (low_link[e.node] > tour_start[node]) {
                    // This is a bridge.
                    is_bridge[e.index] = true;
                    vector<int> component = {node, e.node};
                    sort(component.begin(), component.end());
                    components.push_back(component);
                } else if (low_link[e.node] == tour_start[node]) {
                    // This is the root of a biconnected component.
                    stack.push_back(node);
                    vector<int> component(stack.begin() + size, stack.end());
                    sort(component.begin(), component.end());
                    component.erase(unique(component.begin(), component.end()), component.end());
                    components.push_back(component);
                    stack.resize(size);
                } else {
                    stack.push_back(node);
                }

                // In general, `node` is a cut vertex iff it has a child whose subtree cannot reach above `node`.
                if (low_link[e.node] >= tour_start[node])
                    is_cut[node] = true;
            }
        }

        // The root of the tree is a cut vertex iff it has more than one child.
        if (parent < 0)
            is_cut[node] = children > 1;
    }

    void build(int root = -1) {
        visited.assign(n, false);
        is_cut.assign(n, false);
        is_bridge.assign(edges, false);
        stack.clear();
        components.clear();
        tour = 0;

        if (0 <= root && root < n)
            dfs(root, -1);

        for (int i = 0; i < n; i++)
            if (!visited[i])
                dfs(i, -1);
    }
};

// Note: instead of a block-cut tree this is technically a block-vertex tree, which ends up being much easier to use.
struct block_cut_tree {
    biconnected_components &bi_comps;

    int n, BC, T;
    vector<vector<int>> adj;
    vector<int> parent;
    vector<int> depth;

    // Warning: make sure to call build as well.
    block_cut_tree(biconnected_components &_bi_comps) : bi_comps(_bi_comps) {}

    void dfs(int node, int par) {
        parent[node] = par;
        depth[node] = par < 0 ? 0 : depth[par] + 1;

        for (int neigh : adj[node])
            if (neigh != par)
                dfs(neigh, node);
    }

    void build() {
        n = bi_comps.n;
        BC = int(bi_comps.components.size());
        T = n + BC;
        adj.assign(T, {});

        auto add_edge = [&](int a, int b) {
            assert((a < n) ^ (b < n));
            adj[a].push_back(b);
            adj[b].push_back(a);
        };

        for (int bc = 0; bc < BC; bc++)
            for (int x : bi_comps.components[bc])
                add_edge(x, n + bc);

        parent.assign(T, -1);
        depth.resize(T);

        for (int root = 0; root < T; root++)
            if (parent[root] < 0)
                dfs(root, -1);
    }

    bool same_biconnected_component(int a, int b) const {
        if (depth[a] > depth[b])
            swap(a, b);

        // Two different nodes are in the same biconnected component iff their distance = 2 in the block-cut tree.
        return a == b || (depth[b] == depth[a] + 2 && parent[parent[b]] == a) || (parent[a] >= 0 && parent[a] == parent[b]);
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, M;
    cin >> N >> M;
    biconnected_components bi_comps(N);
    block_cut_tree bc_tree(bi_comps);

    for (int i = 0; i < M; i++) {
        int a, b;
        cin >> a >> b;
        a--; b--;
        bi_comps.add_edge(a, b);
    }

    bi_comps.build();
    bc_tree.build();
    vector<int> cut_vertices;

    for (int i = 0; i < N; i++)
        if (bi_comps.is_cut[i])
            cut_vertices.push_back(i);

    int CV = int(cut_vertices.size());
    cout << CV << '\n';

    for (int i = 0; i < CV; i++)
        cout << cut_vertices[i] + 1 << (i < CV - 1 ? ' ' : '\n');

    for (int i = 0; i < bi_comps.edges; i++)
        if (bi_comps.is_bridge[i])
            cout << bi_comps.edge_list[i][0] + 1 << ' ' << bi_comps.edge_list[i][1] + 1 << '\n';

    int Q;
    cin >> Q;

    for (int q = 0; q < Q; q++) {
        int a, b;
        cin >> a >> b;
        a--; b--;
        cout << bc_tree.same_biconnected_component(a, b);
    }

    cout << '\n';
}
