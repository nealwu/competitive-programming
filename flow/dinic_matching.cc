#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

const int INF = int(1e9) + 5;

struct dinic_matching {
    int n, m;
    vector<vector<int>> adj;
    vector<int> dist;
    vector<int> right_match, left_match;
    vector<int> edge_index;
    bool match_called = false;
    int matches = 0;

    dinic_matching() {
        init(0, 0);
    }

    dinic_matching(int _n, int _m) {
        init(_n, _m);
    }

    void init(int _n, int _m) {
        n = _n;
        m = _m;
        adj.assign(n, {});
        match_called = false;
        matches = 0;
    }

    void add_edge(int a, int b) {
        assert(0 <= a && a < n);
        assert(0 <= b && b < m);
        adj[a].push_back(b);
    }

    bool bfs() {
        vector<int> q(n);
        int q_start = 0, q_end = 0;
        dist.assign(n, INF);

        auto bfs_check = [&](int node, int new_dist) -> void {
            if (new_dist < dist[node]) {
                dist[node] = new_dist;
                q[q_end++] = node;
            }
        };

        for (int i = 0; i < n; i++)
            if (right_match[i] < 0)
                bfs_check(i, 0);

        bool has_path = false;

        while (q_start < q_end) {
            int left = q[q_start++];

            for (int right : adj[left])
                if (left_match[right] < 0)
                    has_path = true;
                else
                    bfs_check(left_match[right], dist[left] + 1);
        }

        return has_path;
    }

    bool dfs(int left) {
        // Because we are only performing DFS in increasing order of dist, we don't have to revisit fully searched edges
        // again later.
        while (edge_index[left] < int(adj[left].size())) {
            int right = adj[left][edge_index[left]++];
            bool solved = false;

            if (left_match[right] < 0 || (dist[left] + 1 == dist[left_match[right]] && dfs(left_match[right]))) {
                left_match[right] = left;
                right_match[left] = right;
                solved = true;
            }

            if (solved)
                return true;
        }

        dist[left] = INF;
        return false;
    }

    int match() {
        match_called = true;
        right_match.assign(n, -1);
        left_match.assign(m, -1);
        matches = 0;

        while (bfs()) {
            edge_index.assign(n, 0);

            for (int i = 0; i < n; i++)
                if (right_match[i] < 0 && dfs(i))
                    matches++;
        }

        return matches;
    }

    vector<bool> reachable_left, reachable_right;

    void _reachable_dfs(int left) {
        reachable_left[left] = true;

        for (int right : adj[left])
            if (right != right_match[left] && !reachable_right[right]) {
                reachable_right[right] = true;
                int next_left = left_match[right];

                if (next_left >= 0 && !reachable_left[next_left])
                    _reachable_dfs(next_left);
            }
    }

    void solve_reachable() {
        reachable_left.assign(n, false);
        reachable_right.assign(m, false);

        for (int i = 0; i < n; i++)
            if (right_match[i] < 0 && !reachable_left[i])
                _reachable_dfs(i);
    }

    // The min vertex cover in a bipartite graph corresponds to the min cut in its flow graph.
    vector<int> min_vertex_cover() {
        assert(match_called);
        solve_reachable();
        vector<int> cover;

        for (int i = 0; i < n; i++)
            if (!reachable_left[i])
                cover.push_back(i);

        for (int i = 0; i < m; i++)
            if (reachable_right[i])
                cover.push_back(n + i);

        assert(int(cover.size()) == matches);
        return cover;
    }

    // The max independent set is the complement of the min vertex cover.
    vector<int> max_independent_set() {
        assert(match_called);
        solve_reachable();
        vector<int> independent_set;

        for (int i = 0; i < n; i++)
            if (reachable_left[i])
                independent_set.push_back(i);

        for (int i = 0; i < m; i++)
            if (!reachable_right[i])
                independent_set.push_back(n + i);

        assert(int(independent_set.size()) + matches == n + m);
        return independent_set;
    }
};


// Accepted on https://www.spoj.com/problems/MATCHING/

int main() {
    int N, M, P;
    scanf("%d %d %d", &N, &M, &P);
    dinic_matching graph(N, M);

    for (int i = 0; i < P; i++) {
        int a, b;
        scanf("%d %d", &a, &b);
        a--; b--;
        graph.add_edge(a, b);
    }

    printf("%d\n", graph.match());
}
