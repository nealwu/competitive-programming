#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

const int INF = int(1e9) + 5;

// Warning: flow_t must be able to handle the sum of flows, not just individual edges.
template<typename flow_t>
struct dinic {
    struct edge {
        int node, rev;
        flow_t capacity, original;

        edge() {}

        edge(int _node, int _rev, flow_t _capacity)
            : node(_node), rev(_rev), capacity(_capacity), original(_capacity) {}
    };

    int V = -1;
    vector<vector<edge>> adj;
    vector<int> dist;
    vector<int> edge_index;
    bool flow_called;

    dinic(int vertices = -1) {
        if (vertices >= 0)
            init(vertices);
    }

    void init(int vertices) {
        V = vertices;
        adj.assign(V, {});
        dist.resize(V);
        edge_index.resize(V);
        flow_called = false;
    }

    void _add_edge(int u, int v, flow_t capacity1, flow_t capacity2) {
        assert(0 <= u && u < V && 0 <= v && v < V);
        assert(capacity1 >= 0 && capacity2 >= 0);
        edge uv_edge(v, int(adj[v].size()) + (u == v ? 1 : 0), capacity1);
        edge vu_edge(u, int(adj[u].size()), capacity2);
        adj[u].push_back(uv_edge);
        adj[v].push_back(vu_edge);
    }

    void add_directional_edge(int u, int v, flow_t capacity) {
        _add_edge(u, v, capacity, 0);
    }

    void add_bidirectional_edge(int u, int v, flow_t capacity) {
        _add_edge(u, v, capacity, capacity);
    }

    edge &reverse_edge(const edge &e) {
        return adj[e.node][e.rev];
    }

    void bfs_check(queue<int> &q, int node, int new_dist) {
        if (new_dist < dist[node]) {
            dist[node] = new_dist;
            q.push(node);
        }
    }

    bool bfs(int source, int sink) {
        dist.assign(V, INF);
        queue<int> q;
        bfs_check(q, source, 0);

        while (!q.empty()) {
            int top = q.front(); q.pop();

            for (edge &e : adj[top])
                if (e.capacity > 0)
                    bfs_check(q, e.node, dist[top] + 1);
        }

        return dist[sink] < INF;
    }

    flow_t dfs(int node, flow_t path_cap, int sink) {
        if (node == sink)
            return path_cap;

        if (dist[node] >= dist[sink])
            return 0;

        flow_t total_flow = 0;

        // Because we are only performing DFS in increasing order of dist, we don't have to revisit fully searched edges
        // again later.
        while (edge_index[node] < int(adj[node].size())) {
            edge &e = adj[node][edge_index[node]];

            if (e.capacity > 0 && dist[node] + 1 == dist[e.node]) {
                flow_t path = dfs(e.node, min(path_cap, e.capacity), sink);
                path_cap -= path;
                e.capacity -= path;
                reverse_edge(e).capacity += path;
                total_flow += path;
            }

            // If path_cap is 0, we don't want to increment edge_index[node] as this edge may not be fully searched yet.
            if (path_cap == 0)
                break;

            edge_index[node]++;
        }

        return total_flow;
    }

    // Can also be used to reverse flow or compute incremental flows after graph modification.
    flow_t flow(int source, int sink, flow_t flow_cap = numeric_limits<flow_t>::max()) {
        assert(V >= 0);
        flow_called = true;
        flow_t total_flow = 0;

        while (flow_cap > 0 && bfs(source, sink)) {
            edge_index.assign(V, 0);
            flow_t increment = dfs(source, flow_cap, sink);
            assert(increment > 0);
            total_flow += increment;
            flow_cap -= increment;
        }

        return total_flow;
    }

    vector<bool> reachable;

    void reachable_dfs(int node) {
        reachable[node] = true;

        for (edge &e : adj[node])
            if (e.capacity > 0 && !reachable[e.node])
                reachable_dfs(e.node);
    }

    void solve_reachable(int source) {
        reachable.assign(V, false);
        reachable_dfs(source);
    }

    // Returns a list of {capacity, {from_node, to_node}} representing edges in the min cut.
    vector<pair<flow_t, pair<int, int>>> min_cut(int source) {
        assert(flow_called);
        solve_reachable(source);
        vector<pair<flow_t, pair<int, int>>> cut;

        for (int node = 0; node < V; node++)
            if (reachable[node])
                for (edge &e : adj[node])
                    if (!reachable[e.node] && e.capacity < e.original)
                        cut.emplace_back(e.original - e.capacity, make_pair(node, e.node));

        return cut;
    }

    // Helper function for setting up incremental / reverse flows. Can become invalid if adding additional edges.
    edge *find_edge(int a, int b) {
        for (edge &e : adj[a])
            if (e.node == b)
                return &e;

        return nullptr;
    }
};


// Accepted on https://www.spoj.com/problems/FASTFLOW/

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, M;
    string str;
    cin >> str;
    bool directed_mode = false;

    if (str == "directed") {
        directed_mode = true;
        cin >> N >> M;
    } else {
        N = stoi(str);
        cin >> M;
    }

    dinic<int64_t> graph(N);

    for (int i = 0; i < M; i++) {
        int a, b, c;
        cin >> a >> b >> c;
        a--; b--;

        if (directed_mode)
            graph.add_directional_edge(a, b, c);
        else
            graph.add_bidirectional_edge(a, b, c);
    }

    int64_t answer = graph.flow(0, N - 1);
    cout << answer << '\n';
    auto cut = graph.min_cut(0);
    int64_t cut_sum = 0;

    for (auto &t : cut)
        cut_sum += t.first;

    assert(answer == cut_sum);
}
