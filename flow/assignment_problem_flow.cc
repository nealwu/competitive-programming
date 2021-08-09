#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

// Warning: when choosing flow_t and cost_t, make sure they can handle the sum of flows and costs, not just individual
// flows and costs.
template<typename flow_t, typename cost_t>
struct min_cost_flow {
    const cost_t COST_INF = numeric_limits<cost_t>::max() / 2;

    struct edge {
        int node, rev;
        flow_t capacity;
        cost_t cost;

        edge() {}

        edge(int _node, int _rev, flow_t _capacity, cost_t _cost)
            : node(_node), rev(_rev), capacity(_capacity), cost(_cost) {}
    };

    int V = -1, E = 0;
    vector<vector<edge>> adj;
    vector<cost_t> dist;
    vector<int> prv;
    vector<edge*> prev_edge;
    bool too_much_bellman_ford = false;

    min_cost_flow(int vertices = -1) {
        if (vertices >= 0)
            init(vertices);
    }

    void init(int vertices) {
        V = vertices;
        E = 0;
        adj.assign(V, {});
        dist.resize(V);
        prv.resize(V);
        prev_edge.resize(V);
        too_much_bellman_ford = false;
    }

    void add_directional_edge(int u, int v, flow_t capacity, cost_t cost) {
        assert(0 <= u && u < V && 0 <= v && v < V);
        assert(capacity >= 0);
        edge uv_edge(v, int(adj[v].size()) + (u == v ? 1 : 0), capacity, cost);
        edge vu_edge(u, int(adj[u].size()), 0, -cost);
        adj[u].push_back(uv_edge);
        adj[v].push_back(vu_edge);
        E++;
    }

    edge &reverse_edge(const edge &e) {
        return adj[e.node][e.rev];
    }

    bool bellman_ford(int source, int sink) {
        dist.assign(V, COST_INF);
        prv.assign(V, -1);
        prev_edge.assign(V, nullptr);

        int64_t work = 0;
        vector<int> last_seen(V, -1);
        vector<int> nodes = {source}, next_nodes;
        dist[source] = 0;

        for (int iteration = 0; iteration < V; iteration++) {
            next_nodes.clear();

            for (int node : nodes) {
                for (edge &e : adj[node])
                    if (e.capacity > 0 && dist[node] + e.cost < dist[e.node]) {
                        dist[e.node] = dist[node] + e.cost;
                        prv[e.node] = node;
                        prev_edge[e.node] = &e;

                        if (last_seen[e.node] != iteration) {
                            last_seen[e.node] = iteration;
                            next_nodes.push_back(e.node);
                        }
                    }

                work += adj[node].size();
            }

            swap(nodes, next_nodes);
        }

        if (work > 1.75L * E * (V == 0 ? 0 : 32 - __builtin_clz(V)) + 100) {
            too_much_bellman_ford = true;
            return false;
        }

        return prv[sink] != -1;
    }

    struct dijkstra_state {
        cost_t dist;
        int node;

        dijkstra_state() {}

        dijkstra_state(cost_t _dist, int _node) : dist(_dist), node(_node) {}

        bool operator<(const dijkstra_state &other) const {
            return dist > other.dist;
        }
    };

    void dijkstra_check(priority_queue<dijkstra_state> &pq, int node, cost_t new_dist, int previous, edge *previous_edge) {
        if (new_dist < dist[node]) {
            dist[node] = new_dist;
            prv[node] = previous;
            prev_edge[node] = previous_edge;
            pq.emplace(dist[node], node);
        }
    }

    bool dijkstra(int source, int sink) {
        dist.assign(V, COST_INF);
        prv.assign(V, -1);
        prev_edge.assign(V, nullptr);

        priority_queue<dijkstra_state> pq;
        dijkstra_check(pq, source, 0, -1, nullptr);

        while (!pq.empty()) {
            dijkstra_state top = pq.top();
            pq.pop();

            if (top.dist > dist[top.node])
                continue;

            for (edge &e : adj[top.node])
                if (e.capacity > 0)
                    dijkstra_check(pq, e.node, top.dist + e.cost, top.node, &e);
        }

        return prv[sink] != -1;
    }

    void reduce_cost() {
        for (int i = 0; i < V; i++)
            for (edge &e : adj[i])
                if (dist[i] < COST_INF && dist[e.node] < COST_INF)
                    e.cost += dist[i] - dist[e.node];
    }

    pair<flow_t, cost_t> solve_min_cost_flow(int source, int sink, flow_t flow_goal = numeric_limits<flow_t>::max()) {
        assert(V >= 0);
        flow_t total_flow = 0;
        cost_t total_cost = 0;
        cost_t reduce_sum = 0;

        auto process_path = [&] {
            flow_t path_cap = flow_goal - total_flow;
            cost_t cost_sum = 0;

            for (int node = sink; prv[node] != -1; node = prv[node])
                path_cap = min(path_cap, prev_edge[node]->capacity);

            for (int node = sink; prv[node] != -1; node = prv[node]) {
                edge *e = prev_edge[node];
                e->capacity -= path_cap;
                reverse_edge(*e).capacity += path_cap;
                cost_sum += e->cost;
            }

            total_flow += path_cap;
            total_cost += (reduce_sum + cost_sum) * path_cap;
        };

        while (total_flow < flow_goal && bellman_ford(source, sink))
            process_path();

        if (too_much_bellman_ford) {
            do {
                reduce_cost();
                reduce_sum += dist[sink];
                process_path();
            } while (total_flow < flow_goal && dijkstra(source, sink));
        }

        return make_pair(total_flow, total_cost);
    }
};

struct assignment_problem {
    int n, m;
    vector<vector<int64_t>> costs;

    assignment_problem(int _n = -1, int _m = -1) : n(_n), m(_m) {
        if (m < 0)
            m = n;

        if (n > 0)
            costs.assign(n, vector<int64_t>(m, 0));
    }

    template<typename T>
    assignment_problem(const vector<vector<T>> &_costs) {
        build(_costs);
    }

    template<typename T>
    void build(const vector<vector<T>> &_costs) {
        n = int(_costs.size());
        m = _costs.empty() ? 0 : int(_costs[0].size());
        costs.assign(n, vector<int64_t>(m, 0));

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                costs[i][j] = _costs[i][j];
    }

    int64_t solve() {
        int v = n + m + 2, source = v - 2, sink = v - 1;
        min_cost_flow<int, int64_t> graph(v);

        for (int i = 0; i < n; i++)
            graph.add_directional_edge(source, i, 1, 0);

        for (int j = 0; j < m; j++)
            graph.add_directional_edge(n + j, sink, 1, 0);

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                graph.add_directional_edge(i, n + j, 1, costs[i][j]);

        return graph.solve_min_cost_flow(source, sink).second;
    }
};


int main() {
    int N, M;
    cin >> N >> M;
    assignment_problem solver;
    vector<vector<int64_t>> costs(N, vector<int64_t>(M, 0));

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
            cin >> costs[i][j];

    solver.build(costs);
    cout << solver.solve() << '\n';
}
