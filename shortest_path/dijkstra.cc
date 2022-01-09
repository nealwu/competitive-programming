#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

const int64_t INF64 = int64_t(2e18) + 5;

struct Dijkstra {
    struct edge {
        int node = -1;
        int64_t weight = 0;

        edge() {}

        edge(int _node, int64_t _weight) : node(_node), weight(_weight) {}
    };

    struct state {
        int64_t dist;
        int node;

        state() {}

        state(int64_t _dist, int _node) : dist(_dist), node(_node) {}

        bool operator<(const state &other) const {
            return dist > other.dist;
        }
    };

    int n;
    vector<vector<edge>> adj;
    vector<int64_t> dist;
    vector<int> parent;

    Dijkstra(int _n = 0) {
        init(_n);
    }

    void init(int _n) {
        n = _n;
        adj.assign(n, {});
    }

    void add_directional_edge(int a, int b, int64_t weight) {
        adj[a].emplace_back(b, weight);
    }

    void add_bidirectional_edge(int a, int b, int64_t weight) {
        add_directional_edge(a, b, weight);
        add_directional_edge(b, a, weight);
    }

    void dijkstra_check(priority_queue<state> &pq, int node, int from, int64_t new_dist) {
        if (new_dist < dist[node]) {
            dist[node] = new_dist;
            parent[node] = from;
            pq.emplace(new_dist, node);
        }
    }

    void dijkstra(const vector<int> &source) {
        if (n == 0) return;
        dist.assign(n, INF64);
        parent.assign(n, -1);
        priority_queue<state> pq;

        for (int src : source)
            dijkstra_check(pq, src, -1, 0);

        while (!pq.empty()) {
            state top = pq.top();
            pq.pop();

            if (top.dist > dist[top.node])
                continue;

            for (edge &e : adj[top.node])
                dijkstra_check(pq, e.node, top.node, top.dist + e.weight);
        }
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, M;
    cin >> N >> M;
    Dijkstra dijkstra(N);

    for (int i = 0; i < M; i++) {
        int a, b;
        int64_t weight;
        cin >> a >> b >> weight;
        a--; b--;
        dijkstra.add_bidirectional_edge(a, b, weight);
    }

    long double begin = clock();
    dijkstra.dijkstra({0});
    fprintf(stderr, "Dijkstra time: %.3Lfs\n", (clock() - begin) / CLOCKS_PER_SEC);

    int64_t total = 0;

    for (int i = 0; i < N; i++)
        total += dijkstra.dist[i];

    cout << total << '\n';
}
