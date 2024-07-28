#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

template<typename A, typename B> ostream& operator<<(ostream &os, const pair<A, B> &p) { return os << '(' << p.first << ", " << p.second << ')'; }
template<typename... Args> ostream& operator<<(ostream& os, const tuple<Args...>& t) { os << '('; apply([&os](const Args&... args) { size_t n = 0; ((os << args << (++n != sizeof...(Args) ? ", " : "")), ...); }, t); return os << ')'; }
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

    struct bi_component {
        vector<int> nodes, edges;
    };

    int n, edges;
    vector<vector<edge>> adj;
    vector<array<int, 2>> edge_list;
    vector<int> tour_start, low_link;
    vector<bool> visited, is_cut;
    vector<bool> is_tree, is_bridge;
    vector<bi_component> components;
    vector<int> edge_to_component;
    vector<int> stack;
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

    int get_deeper_node(int edge_index) const {
        int a = edge_list[edge_index][0], b = edge_list[edge_index][1];
        return tour_start[a] > tour_start[b] ? a : b;
    }

    void create_component(int root, const vector<int> &edge_indices) {
        bi_component component;
        component.edges = edge_indices;

        if (root >= 0)
            component.nodes = {root};

        for (int e_index : edge_indices) {
            if (is_tree[e_index])
                component.nodes.push_back(get_deeper_node(e_index));

            edge_to_component[e_index] = int(components.size());
        }

        components.push_back(component);
    }

    void dfs(int node, int parent_edge) {
        assert(!visited[node]);
        visited[node] = true;
        tour_start[node] = tour++;
        low_link[node] = tour_start[node];
        is_cut[node] = false;
        int children = 0;

        for (edge &e : adj[node]) {
            // Skip the previous edge to the parent, but allow multi-edges.
            if (e.index == parent_edge)
                continue;

            // A self-loop is its own biconnected component.
            if (e.node == node) {
                create_component(-1, {e.index});
                continue;
            }

            if (visited[e.node]) {
                // e.node is a candidate for low_link.
                low_link[node] = min(low_link[node], tour_start[e.node]);

                if (tour_start[e.node] < tour_start[node])
                    stack.push_back(e.index);
            } else {
                // This is a tree edge.
                is_tree[e.index] = true;
                children++;
                size_t size_before = stack.size();
                dfs(e.node, e.index);
                stack.push_back(e.index);

                // e.node is part of our subtree.
                low_link[node] = min(low_link[node], low_link[e.node]);

                if (low_link[e.node] > tour_start[node]) {
                    // This is a bridge.
                    is_bridge[e.index] = true;
                    create_component(node, {e.index});
                    assert(stack.back() == e.index);
                    stack.pop_back();
                } else if (low_link[e.node] == tour_start[node]) {
                    // This is the root of a biconnected component.
                    create_component(node, vector<int>(stack.begin() + size_before, stack.end()));
                    stack.resize(size_before);
                }

                // In general, `node` is a cut vertex iff it has a child whose subtree cannot reach above `node`.
                if (low_link[e.node] >= tour_start[node])
                    is_cut[node] = true;
            }
        }

        // The root of the tree is a cut vertex iff it has more than one child.
        if (parent_edge < 0)
            is_cut[node] = children > 1;
    }

    void build(int root = -1) {
        visited.assign(n, false);
        is_cut.assign(n, false);
        is_tree.assign(edges, false);
        is_bridge.assign(edges, false);
        edge_to_component.assign(edges, -1);
        components.clear();
        stack.clear();
        tour = 0;

        if (0 <= root && root < n)
            dfs(root, -1);

        for (int i = 0; i < n; i++)
            if (!visited[i])
                dfs(i, -1);

        assert(stack.empty());
        assert(find(edge_to_component.begin(), edge_to_component.end(), -1) == edge_to_component.end());
    }
};

// Note: instead of a block-cut tree this is technically a block-vertex tree, which ends up being much easier to use.
struct block_cut_tree {
    biconnected_components &bi_comps;

    int n, B, T;
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
        B = int(bi_comps.components.size());
        T = n + B;
        adj.assign(T, {});

        auto add_edge = [&](int a, int b) -> void {
            assert((a < n) ^ (b < n));
            adj[a].push_back(b);
            adj[b].push_back(a);
        };

        for (int b = 0; b < B; b++)
            for (int x : bi_comps.components[b].nodes)
                add_edge(x, n + b);

        parent.assign(T, -1);
        depth.assign(T, -1);

        for (int root = 0; root < T; root++)
            if (depth[root] < 0)
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
