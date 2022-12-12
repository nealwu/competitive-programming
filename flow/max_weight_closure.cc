#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>
using namespace std;

const int INF = int(1e9) + 5;

enum edge_type : uint8_t { DIRECTIONAL, DIRECTIONAL_REVERSE, BIDIRECTIONAL };

// Warning: flow_t must be able to handle the sum of flows, not just individual edges.
template<typename flow_t>
struct dinic {
    struct edge {
        int node, rev;
        flow_t capacity;
        edge_type type;

        edge() {}

        edge(int _node, int _rev, flow_t _capacity, edge_type _type)
            : node(_node), rev(_rev), capacity(_capacity), type(_type) {}
    };

    int V = -1;
    vector<vector<edge>> adj;
    vector<int> dist;
    vector<int> edge_index;
    bool flow_called = false;

    dinic(int vertices = -1) {
        if (vertices >= 0)
            init(vertices);
    }

    void init(int vertices) {
        V = vertices;
        adj.assign(V, {});
        flow_called = false;
    }

    void _add_edge(int u, int v, flow_t capacity1, flow_t capacity2, edge_type type1, edge_type type2) {
        assert(0 <= min(u, v) && max(u, v) < V);
        assert(capacity1 >= 0 && capacity2 >= 0);
        edge uv_edge(v, int(adj[v].size()) + (u == v ? 1 : 0), capacity1, type1);
        edge vu_edge(u, int(adj[u].size()), capacity2, type2);
        adj[u].push_back(uv_edge);
        adj[v].push_back(vu_edge);
    }

    void add_directional_edge(int u, int v, flow_t capacity) {
        _add_edge(u, v, capacity, 0, DIRECTIONAL, DIRECTIONAL_REVERSE);
    }

    void add_bidirectional_edge(int u, int v, flow_t capacity) {
        _add_edge(u, v, capacity, capacity, BIDIRECTIONAL, BIDIRECTIONAL);
    }

    edge &reverse_edge(const edge &e) {
        return adj[e.node][e.rev];
    }

    bool bfs(int source, int sink) {
        vector<int> q(V);
        int q_start = 0, q_end = 0;
        dist.assign(V, INF);

        auto bfs_check = [&](int node, int new_dist) -> void {
            if (new_dist < dist[node]) {
                dist[node] = new_dist;
                q[q_end++] = node;
            }
        };

        bfs_check(source, 0);

        while (q_start < q_end) {
            int top = q[q_start++];

            for (edge &e : adj[top])
                if (e.capacity > 0)
                    bfs_check(e.node, dist[top] + 1);
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

    void _reachable_dfs(int node) {
        reachable[node] = true;

        for (edge &e : adj[node])
            if (e.capacity > 0 && !reachable[e.node])
                _reachable_dfs(e.node);
    }

    void solve_reachable(int source) {
        reachable.assign(V, false);
        _reachable_dfs(source);
    }

    // Returns a list of {capacity, {from_node, to_node}} representing edges in the min cut.
    vector<pair<flow_t, pair<int, int>>> min_cut(int source) {
        assert(flow_called);
        solve_reachable(source);
        vector<pair<flow_t, pair<int, int>>> cut;

        for (int node = 0; node < V; node++)
            for (edge &e : adj[node])
                if (reachable[node] && !reachable[e.node] && e.type != DIRECTIONAL_REVERSE) {
                    flow_t rev_cap = reverse_edge(e).capacity;
                    flow_t original_cap = e.type == BIDIRECTIONAL ? rev_cap / 2 : rev_cap;
                    cut.emplace_back(original_cap, make_pair(node, e.node));
                }

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

// max_weight_closure solves the following problem:
// There are n projects you can complete. The i-th project gives you a gain (or loss if negative) of P_i money.
// Projects can have dependencies on other projects; you can do a project as long as you also do all its dependencies.
// Cyclic dependencies are allowed; if a group of projects is strongly connected, you have to do all or none of them.
// See https://en.wikipedia.org/wiki/Closure_problem. "Project A depends on project B" means an A -> B edge.
// Warning: cost_t must be able to handle the sum of costs, not just individual amounts.
template<typename cost_t>
struct max_weight_closure {
    int n, source, sink;
    dinic<cost_t> graph;
    vector<bool> chosen;
    cost_t positive_total = 0;

    max_weight_closure() {}

    max_weight_closure(int _n) {
        init(_n);
    }

    template<typename T_array>
    max_weight_closure(const T_array &projects) {
        init(projects);
    }

    void init(int _n) {
        n = _n;
        int V = n + 2;
        source = V - 2;
        sink = V - 1;
        graph.init(V);
    }

    template<typename T_array>
    void init(const T_array &projects) {
        init(int(projects.size()));
        set_projects(projects);
    }

    template<typename T_array>
    void set_projects(const T_array &projects) {
        chosen.assign(n, false);
        positive_total = 0;

        for (int i = 0; i < n; i++)
            if (projects[i] >= 0) {
                graph.add_directional_edge(source, i, projects[i]);
                positive_total += projects[i];
                chosen[i] = true;
            } else {
                graph.add_directional_edge(i, sink, -projects[i]);
            }
    }

    // Project `a` depends on project `b`.
    void add_dependency(int a, int b) {
        assert(0 <= min(a, b) && max(a, b) < n);
        graph.add_directional_edge(a, b, numeric_limits<cost_t>::max());
    }

    cost_t solve() {
        return positive_total - graph.flow(source, sink);
    }

    vector<int> chosen_projects() {
        auto cut = graph.min_cut(source);

        for (auto &cut_edge : cut)
            if (cut_edge.second.first == source)
                chosen[cut_edge.second.second] = false;
            else if (cut_edge.second.second == sink)
                chosen[cut_edge.second.first] = true;

        vector<int> solution;

        for (int i = 0; i < n; i++)
            if (chosen[i])
                solution.push_back(i);

        return solution;
    }
};


template<typename A, typename B> ostream& operator<<(ostream &os, const pair<A, B> &p) { return os << '(' << p.first << ", " << p.second << ')'; }
template<typename T_container, typename T = typename enable_if<!is_same<T_container, string>::value, typename T_container::value_type>::type> ostream& operator<<(ostream &os, const T_container &v) { os << '{'; string sep; for (const T &x : v) os << sep << x, sep = ", "; return os << '}'; }

void dbg_out() { cerr << endl; }
template<typename Head, typename... Tail> void dbg_out(Head H, Tail... T) { cerr << ' ' << H; dbg_out(T...); }
#ifdef NEAL_DEBUG
#define dbg(...) cerr << '[' << __FILE__ << ':' << __LINE__ << "] (" << #__VA_ARGS__ << "):", dbg_out(__VA_ARGS__)
#else
#define dbg(...)
#endif

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, M;
    cin >> N >> M;
    max_weight_closure<int64_t> solver(N);
    vector<int> projects(N);

    for (auto &p : projects)
        cin >> p;

    solver.set_projects(projects);
    vector<array<int, 2>> edges(M);

    for (int i = 0; i < M; i++) {
        int a, b;
        cin >> a >> b;
        a--; b--;
        edges[i] = {a, b};
        solver.add_dependency(a, b);
    }

    int64_t ans = solver.solve();
    cout << ans << '\n';

    vector<int> solution = solver.chosen_projects();
    int64_t sum = 0;

    for (int p : solution)
        sum += projects[p];

    assert(ans == sum);

    for (auto &e : edges)
        assert(!(solver.chosen[e[0]] && !solver.chosen[e[1]]));
}
