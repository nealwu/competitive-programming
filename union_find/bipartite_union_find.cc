#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;

struct bipartite_union_find {
    // When data[x] < 0, x is a root and -data[x] is its tree size. When data[x] >= 0, data[x] is x's parent.
    vector<int> data;
    vector<bool> bipartite;  // whether the component rooted at x is bipartite
    vector<bool> edge_parity;  // the parity of the edge x -- parent[x]
    int components = 0;

    bipartite_union_find(int n = -1) {
        if (n >= 0)
            init(n);
    }

    void init(int n) {
        data.assign(n, -1);
        bipartite.assign(n, true);
        edge_parity.assign(n, false);
        components = n;
    }

    int find(int x) {
        if (data[x] < 0)
            return x;

        int root = find(data[x]);
        edge_parity[x] = edge_parity[x] ^ edge_parity[data[x]];
        return data[x] = root;
    }

    int get_size(int x) {
        return -data[find(x)];
    }

    // Returns true if x and y are in the same component.
    bool query_component(int x, int y) {
        return find(x) == find(y);
    }

    // Returns the parity status between x and y (0 = same, 1 = different). Requires them to be in the same component.
    bool query_parity(int x, int y) {
        bool same_component = query_component(x, y);
        assert(same_component);
        return edge_parity[x] ^ edge_parity[y];
    }

    // Returns {union succeeded, edge consistent with bipartite conditions}.
    pair<bool, bool> unite(int x, int y, bool different = true) {
        int x_root = find(x);
        int y_root = find(y);
        bool root_parity = edge_parity[x] ^ edge_parity[y] ^ different;
        x = x_root;
        y = y_root;

        if (x == y) {
            bool consistent = !root_parity;
            bipartite[x] = bipartite[x] && consistent;
            return {false, consistent};
        }

        if (-data[x] < -data[y])
            swap(x, y);

        data[x] += data[y];
        data[y] = x;
        bipartite[x] = bipartite[x] && bipartite[y];
        edge_parity[y] = root_parity;
        components--;
        return {true, true};
    }

    // Add an assertion that x and y are different; i.e., a normal edge.
    pair<bool, bool> add_different_edge(int x, int y) {
        return unite(x, y, true);
    }

    // Add an assertion that x and y are the same.
    pair<bool, bool> add_same_edge(int x, int y) {
        return unite(x, y, false);
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, Q;
    cin >> N >> Q;
    bipartite_union_find UF(N);

    for (int q = 0; q < Q; q++) {
        int type, a, b;
        cin >> type >> a >> b;
        assert(1 <= min(a, b) && max(a, b) <= N);
        a--; b--;

        if (type == 1) {
            bool same_component = UF.query_component(a, b);

            if (same_component)
                cout << same_component << ' ' << UF.query_parity(a, b) << '\n';
            else
                cout << same_component << '\n';
        } else if (type == 2) {
            int e;
            cin >> e;
            pair<bool, bool> result = UF.unite(a, b, e);
            cout << result.first << ' ' << result.second << '\n';
        } else {
            assert(false);
        }
    }

    for (int i = 0; i < N; i++)
        cout << UF.bipartite[UF.find(i)];

    cout << '\n';
}
