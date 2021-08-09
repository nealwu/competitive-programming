// Solution to https://leetcode.com/contest/weekly-contest-154/problems/critical-connections-in-a-network/
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

template<typename T> ostream& operator<<(ostream &os, const vector<T> &v) { os << "{"; string sep; for (const auto &x : v) os << sep << x, sep = ", "; return os << "}"; }

struct find_bridges {
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
    vector<bool> is_bridge;
    int tour;

    find_bridges(int _n = 0) {
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
        int parent_count = 0;

        for (edge &e : adj[node]) {
            // Skip the first edge to the parent, but allow multi-edges.
            if (e.node == parent && parent_count++ == 0)
                continue;

            if (visited[e.node]) {
                // e.node is a candidate for low_link.
                low_link[node] = min(low_link[node], tour_start[e.node]);
            } else {
                dfs(e.node, node);
                is_bridge[e.index] = low_link[e.node] > tour_start[node];
                // e.node is part of our subtree.
                low_link[node] = min(low_link[node], low_link[e.node]);
            }
        }
    }

    void solve() {
        visited.assign(n, false);
        is_bridge.assign(edges, false);
        tour = 0;

        for (int i = 0; i < n; i++)
            if (!visited[i])
                dfs(i, -1);
    }
};


class Solution {
public:
    vector<vector<int>> criticalConnections(int n, vector<vector<int>>& connections) {
        find_bridges graph(n);
        int E = int(connections.size());

        for (int e = 0; e < E; e++) {
            int a = connections[e][0], b = connections[e][1];
            graph.add_edge(a, b);
        }

        graph.solve();

        // Construct the bridge list in the same order as given in the input.
        vector<vector<int>> bridge_list;

        for (int e = 0; e < E; e++)
            if (graph.is_bridge[e])
                bridge_list.push_back(connections[e]);

        return bridge_list;
    }
};

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    vector<vector<int>> connections;
    vector<vector<int>> desired;

    connections = {{0, 1}, {1, 2}, {2, 0}, {1, 3}};
    desired = {{1, 3}};
    assert(Solution().criticalConnections(4, connections) == desired);

    connections = {{0, 1}, {1, 2}, {1, 3}, {4, 3}};
    desired = {{0, 1}, {1, 2}, {1, 3}, {4, 3}};
    assert(Solution().criticalConnections(5, connections) == desired);

    connections = {{0, 1}, {1, 2}, {2, 0}, {1, 3}, {3, 4}, {4, 5}, {5, 3}};
    desired = {{1, 3}};
    assert(Solution().criticalConnections(6, connections) == desired);

    cerr << "Tests passed!" << endl;

    int n, m;
    cin >> n >> m;
    connections = {};

    for (int i = 0; i < m; i++) {
        int a, b;
        cin >> a >> b;
        connections.push_back({a, b});
    }

    vector<vector<int>> answer = Solution().criticalConnections(n, connections);

    for (auto &e : answer) {
        assert(e.size() == 2);
        cout << e[0] << ' ' << e[1] << '\n';
    }
}
