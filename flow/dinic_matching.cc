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
    vector<bool> matched;
    vector<int> edge_index;
    vector<int> prv;

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
        dist.resize(n);
        matched.resize(n);
        edge_index.resize(n);
        prv.resize(m);
    }

    void add_edge(int a, int b) {
        assert(0 <= a && a < n);
        assert(0 <= b && b < m);
        adj[a].push_back(b);
    }

    void bfs_check(queue<int> &q, int left, int new_dist) {
        if (new_dist < dist[left]) {
            dist[left] = new_dist;
            q.push(left);
        }
    }

    bool bfs() {
        dist.assign(n, INF);
        queue<int> q;

        for (int i = 0; i < n; i++)
            if (!matched[i])
                bfs_check(q, i, 0);

        bool has_path = false;

        while (!q.empty()) {
            int left = q.front(); q.pop();

            for (int right : adj[left])
                if (prv[right] < 0)
                    has_path = true;
                else
                    bfs_check(q, prv[right], dist[left] + 1);
        }

        return has_path;
    }

    bool dfs(int left) {
        // Because we are only performing DFS in increasing order of dist, we don't have to revisit fully searched edges
        // again later.
        while (edge_index[left] < int(adj[left].size())) {
            int right = adj[left][edge_index[left]++];
            bool solved = false;

            if (prv[right] < 0 || (dist[left] + 1 == dist[prv[right]] && dfs(prv[right]))) {
                prv[right] = left;
                matched[left] = true;
                solved = true;
            }

            if (solved)
                return true;
        }

        dist[left] = INF;
        return false;
    }

    int match() {
        matched.assign(n, false);
        prv.assign(m, -1);
        int matches = 0;

        while (bfs()) {
            edge_index.assign(n, 0);

            for (int i = 0; i < n; i++)
                if (!matched[i] && dfs(i))
                    matches++;
        }

        return matches;
    }

    vector<bool> reachable_left, reachable_right;

    void reachable_dfs(int left) {
        reachable_left[left] = true;

        for (int right : adj[left])
            if (prv[right] != left && !reachable_right[right]) {
                reachable_right[right] = true;
                int next_left = prv[right];

                if (next_left >= 0 && !reachable_left[next_left])
                    reachable_dfs(next_left);
            }
    }

    // The min vertex cover in a bipartite graph corresponds to the min cut in its flow graph.
    vector<int> min_vertex_cover() {
        int match_size = match();
        reachable_left.assign(n, false);
        reachable_right.assign(m, false);

        for (int i = 0; i < n; i++)
            if (!matched[i] && !reachable_left[i])
                reachable_dfs(i);

        vector<int> cover;

        for (int i = 0; i < n; i++)
            if (!reachable_left[i])
                cover.push_back(i);

        for (int i = 0; i < m; i++)
            if (reachable_right[i])
                cover.push_back(n + i);

        assert(int(cover.size()) == match_size);
        return cover;
    }

    // The max independent set is the complement of the min vertex cover.
    vector<int> max_independent_set() {
        int cover_size = int(min_vertex_cover().size());
        vector<int> independent_set;

        for (int i = 0; i < n; i++)
            if (reachable_left[i])
                independent_set.push_back(i);

        for (int i = 0; i < m; i++)
            if (!reachable_right[i])
                independent_set.push_back(n + i);

        assert(int(independent_set.size()) + cover_size == n + m);
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
