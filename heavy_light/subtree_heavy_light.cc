#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

struct segment_change {
    // Use a sentinel value rather than a boolean to save significant memory (4-8 bytes per object).
    static const int SENTINEL = numeric_limits<int>::lowest();

    // Note that to_set goes first, and to_add goes after.
    // TODO: check if these values can overflow int.
    int to_set, to_add;

    // TODO: make sure the default constructor is the identity segment_change.
    segment_change(int _to_add = 0, int _to_set = SENTINEL) : to_set(_to_set), to_add(_to_add) {}

    bool has_set() const {
        return to_set != SENTINEL;
    }

    bool has_change() const {
        return has_set() || to_add != 0;
    }

    // Return the combined result of applying this segment_change followed by `other`.
    // TODO: make sure to check for sentinel values.
    segment_change combine(const segment_change &other) const {
        if (other.has_set())
            return other;

        return segment_change(to_add + other.to_add, to_set);
    }
};

struct segment {
    // TODO: check if these values can overflow int.
    int maximum;
    int64_t sum;

    // TODO: make sure the default constructor is the identity segment.
    segment(int _maximum = numeric_limits<int>::lowest(), int64_t _sum = 0) : maximum(_maximum), sum(_sum) {}

    bool empty() const {
        return maximum == numeric_limits<int>::lowest();
    }

    void apply(int length, const segment_change &change) {
        if (change.has_set()) {
            maximum = change.to_set;
            sum = int64_t(length) * change.to_set;
        }

        maximum += change.to_add;
        sum += int64_t(length) * change.to_add;
    }

    void join(const segment &other) {
        if (empty()) {
            *this = other;
            return;
        } else if (other.empty()) {
            return;
        }

        maximum = max(maximum, other.maximum);
        sum += other.sum;
    }

    // TODO: decide whether to re-implement this for better performance. Mainly relevant when segments contain arrays.
    void join(const segment &a, const segment &b) {
        *this = a;
        join(b);
    }
};

pair<int, int> right_half[32];

struct seg_tree {
    static int highest_bit(unsigned x) {
        return x == 0 ? -1 : 31 - __builtin_clz(x);
    }

    int tree_n = 0;
    vector<segment> tree;
    vector<segment_change> changes;

    seg_tree(int n = -1) {
        if (n >= 0)
            init(n);
    }

    void init(int n) {
        tree_n = 1;

        while (tree_n < n)
            tree_n *= 2;

        tree.assign(2 * tree_n, segment());
        changes.assign(tree_n, segment_change());
    }

    // Builds our tree from an array in O(n).
    void build(const vector<segment> &initial) {
        int n = int(initial.size());
        init(n);
        assert(n <= tree_n);

        for (int i = 0; i < n; i++)
            tree[tree_n + i] = initial[i];

        for (int position = tree_n - 1; position > 0; position--)
            tree[position].join(tree[2 * position], tree[2 * position + 1]);
    }

    void apply_and_combine(int position, int length, const segment_change &change) {
        tree[position].apply(length, change);

        if (position < tree_n)
            changes[position] = changes[position].combine(change);
    }

    void push_down(int position, int length) {
        if (changes[position].has_change()) {
            apply_and_combine(2 * position, length / 2, changes[position]);
            apply_and_combine(2 * position + 1, length / 2, changes[position]);
            changes[position] = segment_change();
        }
    }

    // Calls push_down for all necessary nodes in order to query the range [a, b).
    void push_all(int a, int b) {
        assert(0 <= a && a < b && b <= tree_n);
        a += tree_n;
        b += tree_n - 1;

        for (int up = highest_bit(tree_n); up > 0; up--) {
            int x = a >> up, y = b >> up;
            push_down(x, 1 << up);
            if (x != y) push_down(y, 1 << up);
        }
    }

    void join_and_apply(int position, int length) {
        tree[position].join(tree[2 * position], tree[2 * position + 1]);
        tree[position].apply(length, changes[position]);
    }

    // Calls join for all necessary nodes after updating the range [a, b).
    void join_all(int a, int b) {
        assert(0 <= a && a < b && b <= tree_n);
        a += tree_n;
        b += tree_n - 1;
        int length = 1;

        while (a > 1) {
            a /= 2;
            b /= 2;
            length *= 2;
            join_and_apply(a, length);
            if (a != b) join_and_apply(b, length);
        }
    }

    template<typename T_range_op>
    void process_range(int a, int b, bool needs_join, T_range_op &&range_op) {
        if (a == b) return;
        push_all(a, b);
        int original_a = a, original_b = b;
        int length = 1, r_size = 0;

        for (a += tree_n, b += tree_n; a < b; a /= 2, b /= 2, length *= 2) {
            if (a & 1)
                range_op(a++, length);

            if (b & 1)
                right_half[r_size++] = {--b, length};
        }

        for (int i = r_size - 1; i >= 0; i--)
            range_op(right_half[i].first, right_half[i].second);

        if (needs_join)
            join_all(original_a, original_b);
    }

    segment query(int a, int b) {
        assert(0 <= a && a <= b && b <= tree_n);
        segment answer;

        process_range(a, b, false, [&](int position, int) {
            answer.join(tree[position]);
        });

        return answer;
    }

    segment query_full() const {
        return tree[1];
    }

    void update(int a, int b, const segment_change &change) {
        assert(0 <= a && a <= b && b <= tree_n);

        process_range(a, b, true, [&](int position, int length) {
            apply_and_combine(position, length, change);
        });
    }
};


struct subtree_heavy_light {
    int n = 0;
    bool vertex_mode;

    vector<vector<int>> adj;
    vector<int> parent, depth, subtree_size;

    vector<int> tour_start, tour_end;
    vector<int> chain_root;
    seg_tree full_tree;

    subtree_heavy_light() {}

    subtree_heavy_light(int _n, bool _vertex_mode) {
        init(_n, _vertex_mode);
    }

    void init(int _n, bool _vertex_mode) {
        n = _n;
        vertex_mode = _vertex_mode;

        adj.assign(n, {});
        parent.resize(n);
        depth.resize(n);
        subtree_size.resize(n);

        tour_start.resize(n);
        tour_end.resize(n);
        chain_root.resize(n);
    }

    void add_edge(int a, int b) {
         adj[a].push_back(b);
         adj[b].push_back(a);
    }

    void dfs(int node, int par) {
        parent[node] = par;
        depth[node] = par < 0 ? 0 : depth[par] + 1;
        subtree_size[node] = 1;

        // Erase the edge to parent.
        adj[node].erase(remove(adj[node].begin(), adj[node].end(), par), adj[node].end());

        for (int child : adj[node]) {
            dfs(child, node);
            subtree_size[node] += subtree_size[child];
        }

        // Heavy-light subtree reordering.
        sort(adj[node].begin(), adj[node].end(), [&](int a, int b) {
            return subtree_size[a] > subtree_size[b];
        });
    }

    int tour;

    void chain_dfs(int node, bool heavy) {
        chain_root[node] = heavy ? chain_root[parent[node]] : node;
        tour_start[node] = tour++;
        bool heavy_child = true;

        for (int child : adj[node]) {
            chain_dfs(child, heavy_child);
            heavy_child = false;
        }

        tour_end[node] = tour;
    }

    void build(const segment &initial) {
        tour = 0;
        parent.assign(n, -1);

        for (int i = 0; i < n; i++)
            if (parent[i] < 0) {
                dfs(i, -1);
                chain_dfs(i, false);
            }

        full_tree.init(n);
        full_tree.build(vector<segment>(n, initial));
    }

    void update_subtree(int v, const segment_change &change) {
        full_tree.update(tour_start[v] + (vertex_mode ? 0 : 1), tour_end[v], change);
    }

    segment query_subtree(int v) {
        return full_tree.query(tour_start[v] + (vertex_mode ? 0 : 1), tour_end[v]);
    }

    template<typename T_tree_op>
    int process_path(int u, int v, T_tree_op &&op) {
        while (chain_root[u] != chain_root[v]) {
            // Always pull up the chain with the deeper root.
            if (depth[chain_root[u]] > depth[chain_root[v]])
                swap(u, v);

            int root = chain_root[v];
            op(full_tree, tour_start[root], tour_start[v] + 1);
            v = parent[root];
        }

        if (depth[u] > depth[v])
            swap(u, v);

        // u is now an ancestor of v.
        op(full_tree, tour_start[u] + (vertex_mode ? 0 : 1), tour_start[v] + 1);
        return u;
    }

    int get_lca(int u, int v) {
        return process_path(u, v, [&](seg_tree &, int, int){});
    }

    segment query_path(int u, int v) {
        segment answer;

        process_path(u, v, [&](seg_tree &tree, int a, int b) {
            answer.join(tree.query(a, b));
        });

        return answer;
    }

    void update_path(int u, int v, const segment_change &change) {
        process_path(u, v, [&](seg_tree &tree, int a, int b) {
            tree.update(a, b, change);
        });
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, Q;
    bool vertex_mode;
    cin >> N >> Q >> vertex_mode;
    subtree_heavy_light HLD(N, vertex_mode);
    vector<pair<int, int>> edges;

    for (int i = 0; i < N - 1; i++) {
        int a, b;
        cin >> a >> b;
        a--; b--;
        HLD.add_edge(a, b);
        edges.emplace_back(a, b);
    }

    HLD.build(segment(0, 0));

    for (int q = 0; q < Q; q++) {
        int type, a, b, x;
        cin >> type >> a;
        a--;

        if (type <= 4) {
            cin >> b;
            b--;
        }

        if (type % 4 == 1 || type % 4 == 2)
            cin >> x;

        // There are eight types of operations:
        // 1) Add x to all edges (vertices) on the path from a to b.
        // 2) Set all edges (vertices) on the path from a to b to x.
        // 3) Find the maximum of all edges (vertices) on the path from a to b.
        // 4) Find the sum of all edges (vertices) on the path from a to b.
        // 5) Add x to all edges (vertices) in the subtree rooted at a.
        // 6) Set all edges (vertices) in the subtree rooted at a to x.
        // 7) Find the maximum of all edges (vertices) in the subtree rooted at a.
        // 8) Find the sum of all edges (vertices) in the subtree rooted at a.

        if (type == 1) {
            HLD.update_path(a, b, segment_change(x));
        } else if (type == 2) {
            HLD.update_path(a, b, segment_change(0, x));
        } else if (type <= 4) {
            assert(type == 3 || type == 4);
            segment path = HLD.query_path(a, b);
            cout << (type == 3 ? max(path.maximum, -1) : path.sum) << '\n';
        } else if (type == 5) {
            HLD.update_subtree(a, segment_change(x));
        } else if (type == 6) {
            HLD.update_subtree(a, segment_change(0, x));
        } else {
            assert(type == 7 || type == 8);
            segment subtree = HLD.query_subtree(a);
            cout << (type == 7 ? max(subtree.maximum, -1) : subtree.sum) << '\n';
        }
    }
}
