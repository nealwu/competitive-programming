#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

// A persistent array that uses log n time per lookup, as well as log n time and log n memory per update.
// TODO: persistent trees have both bad constant factor and high memory usage. Make sure we won't TLE or MLE.
template<typename T>
struct persistent_array {
    int tree_n = 0;
    vector<array<int, 2>> tree;
    vector<T> values;
    int tree_reserve_size = INT32_MAX, value_reserve_size = INT32_MAX;

    persistent_array(int n = -1, int max_updates = 0) {
        if (n >= 0)
            init(n, max_updates);
    }

    persistent_array(const vector<T> &v, int max_updates = 0) {
        init(v, max_updates);
    }

    void init(int n, int max_updates = 0) {
        init(vector<T>(n, 0), max_updates);
    }

    // Directly assigning `tree[position] = f(...)` results in segmentation faults, because the address for
    // `tree[position]` can be computed before calling `f()`, which may reallocate `tree`.
    void set_children(int position, int left, int right) { tree[position] = {left, right}; }
    void set_left(int position, int left) { tree[position][0] = left; }
    void set_right(int position, int right) { tree[position][1] = right; }

    int build_tree(int start, int end) {
        if (start >= end)
            return -1;

        int node = int(tree.size());
        tree.emplace_back();

        // At leaves, we set `tree[node]` to point to the appropriate index in `values`.
        if (end - start == 1) {
            set_children(node, start, start);
            return node;
        }

        int mid = (start + end) / 2;
        set_children(node, build_tree(start, mid), build_tree(mid, end));
        return node;
    }

    void init(const vector<T> &v, int max_updates = 0) {
        tree_n = int(v.size());
        values = v;
        tree = {{-1, -1}};
        tree_reserve_size = value_reserve_size = INT32_MAX;

        if (max_updates > 0) {
            // We need to add one to tree_height if tree_n is not a power of two.
            int tree_height = 32 - __builtin_clz(tree_n) + ((tree_n & (tree_n - 1)) != 0);
            tree_reserve_size = 2 * tree_n + max_updates * tree_height;
            value_reserve_size = tree_n + max_updates;
            tree.reserve(tree_reserve_size);
            values.reserve(value_reserve_size);
        }

        build_tree(0, tree_n);
    }

    T get(int root, int index) const {
        assert(root > 0 && 0 <= index && index < tree_n);
        int current = root;
        int start = 0, end = tree_n;

        while (end - start > 1) {
            int mid = (start + end) / 2;

            if (index < mid) {
                current = tree[current][0];
                end = mid;
            } else {
                current = tree[current][1];
                start = mid;
            }
        }

        assert(start == index && end == index + 1);
        return values[tree[current][0]];
    }

    int make_copy(int position) {
        assert(0 <= position && position < int(tree.size()));
        tree.push_back(tree[position]);
        assert(int(tree.size()) <= tree_reserve_size);
        return int(tree.size()) - 1;
    }

    int update_tree(int position, int start, int end, int index, T value) {
        assert(start < end);
        position = make_copy(position);

        if (end - start == 1) {
            assert(start == index);
            set_children(position, int(values.size()), int(values.size()));
            values.push_back(value);
            assert(int(values.size()) <= value_reserve_size);
            return position;
        }

        int mid = (start + end) / 2;

        if (index < mid)
            set_left(position, update_tree(tree[position][0], start, mid, index, value));
        else
            set_right(position, update_tree(tree[position][1], mid, end, index, value));

        return position;
    }

    int update(int root, int index, T value) {
        assert(root > 0 && 0 <= index && index < tree_n);
        return update_tree(root, 0, tree_n, index, value);
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, Q;
    cin >> N >> Q;
    persistent_array<int> values(N, Q);
    vector<int> roots(Q + 1, 1);

    for (int q = 1; q <= Q; q++) {
        string type;
        int index, value;
        cin >> type >> index;

        if (type == "load") {
            roots[q] = roots[index];
        } else if (type == "set") {
            index--;
            cin >> value;
            roots[q] = values.update(roots[q - 1], index, value);
        } else if (type == "query") {
            index--;
            roots[q] = roots[q - 1];
            cout << values.get(roots[q], index) << '\n';
        }
    }
}
