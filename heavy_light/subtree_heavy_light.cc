#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <limits>
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
    void join(const segment &seg0, const segment &seg1) {
        *this = seg0;
        join(seg1);
    }
};

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

        tree.assign(2 * tree_n, {});
        changes.assign(tree_n, {});
    }

    // Builds our tree from an array in O(n).
    void build(const vector<segment> &initial) {
        int n = int(initial.size());
        init(n);
        copy(initial.begin(), initial.end(), tree.begin() + tree_n);

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
        int last = max(highest_bit(a ^ b), 0);

        for (int up = highest_bit(tree_n); up > last; up--)
            push_down(a >> up, 1 << up);

        for (int up = last; up > 0; up--) {
            push_down(a >> up, 1 << up);
            push_down(b >> up, 1 << up);
        }
    }

    void join_and_apply(int position, int length) {
        tree[position].join(tree[2 * position], tree[2 * position + 1]);
        tree[position].apply(length, changes[position]);
    }

    // Calls join_and_apply for all necessary nodes after updating the range [a, b).
    void join_all(int a, int b) {
        assert(0 <= a && a < b && b <= tree_n);
        a += tree_n;
        b += tree_n - 1;
        int last = max(highest_bit(a ^ b), 0);

        for (int up = 1; up <= last; up++) {
            join_and_apply(a >> up, 1 << up);
            join_and_apply(b >> up, 1 << up);
        }

        for (int up = last + 1; up <= highest_bit(tree_n); up++)
            join_and_apply(a >> up, 1 << up);
    }

    template<typename T_range_op>
    void process_range(int a, int b, bool needs_join, T_range_op &&range_op) {
        assert(0 <= a && a <= b && b <= tree_n);
        if (a == b) return;
        push_all(a, b);
        int original_a = a, original_b = b;
        a += tree_n;
        b += tree_n;
        a--;
        int anc_depth = highest_bit(a ^ b);
        int anc_mask = (1 << anc_depth) - 1;

        // Iterate the 0-bits of `a` bottom-up and the 1-bits of `b` top-down.
        for (int v = ~a & anc_mask; v != 0; v &= v - 1) {
            int i = __builtin_ctz(v);
            if (range_op((a >> i) + 1, 1 << i)) return;  // Early return used for search only
        }

        for (int v = b & anc_mask, i; v != 0; v ^= 1 << i) {
            i = highest_bit(v);
            if (range_op((b >> i) - 1, 1 << i)) return;
        }

        if (needs_join)
            join_all(original_a, original_b);
    }

    segment query(int a, int b) {
        segment answer;

        process_range(a, b, false, [&](int position, int) -> bool {
            answer.join(tree[position]);
            return false;
        });

        return answer;
    }

    segment query_full() const {
        return tree[1];
    }

    segment query_single(int index) {
        assert(0 <= index && index < tree_n);
        int position = tree_n + index;

        for (int up = highest_bit(tree_n); up > 0; up--)
            push_down(position >> up, 1 << up);

        return tree[position];
    }

    void update(int a, int b, const segment_change &change) {
        process_range(a, b, true, [&](int position, int length) -> bool {
            apply_and_combine(position, length, change);
            return false;
        });
    }

    void update_single(int index, const segment &seg) {
        assert(0 <= index && index < tree_n);
        int position = tree_n + index;

        for (int up = highest_bit(tree_n); up > 0; up--)
            push_down(position >> up, 1 << up);

        tree[position] = seg;

        while (position > 1) {
            position /= 2;
            tree[position].join(tree[2 * position], tree[2 * position + 1]);
        }
    }

    vector<segment> to_array(int n) {
        for (int i = 1; i < tree_n; i++)
            push_down(i, tree_n >> highest_bit(i));

        return vector<segment>(tree.begin() + tree_n, tree.begin() + tree_n + n);
    }

    // Finds the end of the last subarray starting at `first` satisfying `should_join` via binary search in O(log n).
    template<typename T_bool>
    int find_last_subarray(T_bool &&should_join, int n, int first = 0) {
        assert(0 <= first && first <= n);
        segment current;

        // Check the degenerate case.
        if (!should_join(current, current))
            return first - 1;

        int node = -1;

        // Try to build the range [first, n); when a node fails, search down instead.
        process_range(first, n, false, [&](int position, int) -> bool {
            if (should_join(current, tree[position])) {
                current.join(tree[position]);
                return false;
            }

            node = position;
            return true;
        });

        if (node < 0)
            return n;

        while (node < tree_n) {
            push_down(node, tree_n >> highest_bit(node));

            if (should_join(current, tree[2 * node])) {
                current.join(tree[2 * node]);
                node = 2 * node + 1;
            } else {
                node = 2 * node;
            }
        }

        return node - tree_n;
    }
};


struct subtree_heavy_light {
    int n = 0;
    bool vertex_mode = false;

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
        sort(adj[node].begin(), adj[node].end(), [&](int a, int b) -> bool {
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

    segment query_subtree(int v) {
        return full_tree.query(tour_start[v] + (vertex_mode ? 0 : 1), tour_end[v]);
    }

    void update_subtree(int v, const segment_change &change) {
        full_tree.update(tour_start[v] + (vertex_mode ? 0 : 1), tour_end[v], change);
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
        return process_path(u, v, [&](seg_tree &, int, int) -> void {});
    }

    segment query_path(int u, int v) {
        segment answer;

        process_path(u, v, [&](seg_tree &tree, int a, int b) -> void {
            answer.join(tree.query(a, b));
        });

        return answer;
    }

    void update_path(int u, int v, const segment_change &change) {
        process_path(u, v, [&](seg_tree &tree, int a, int b) -> void {
            tree.update(a, b, change);
        });
    }

    void update_single(int v, const segment &seg) {
        // If vertex mode we update the node; otherwise we update the edge up from the node.
        full_tree.update_single(tour_start[v], seg);
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, Q;
    bool vertex_mode = false;
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
            if ((vertex_mode && a == b) || (!vertex_mode && HLD.parent[a] == b))
                HLD.update_single(a, segment(x, x));
            else
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
