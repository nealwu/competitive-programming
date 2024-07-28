#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>
using namespace std;

// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2016/p0200r0.html
template<class Fun> class y_combinator_result {
    Fun fun_;
public:
    template<class T> explicit y_combinator_result(T &&fun): fun_(std::forward<T>(fun)) {}
    template<class ...Args> decltype(auto) operator()(Args &&...args) { return fun_(std::ref(*this), std::forward<Args>(args)...); }
};
template<class Fun> decltype(auto) y_combinator(Fun &&fun) { return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun)); }


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
    int first, last, max_diff;

    // TODO: make sure the default constructor is the identity segment.
    segment(int _maximum = numeric_limits<int>::lowest(), int64_t _sum = 0, int _first = 0, int _last = 0,
            int _max_diff = -1) : maximum(_maximum), sum(_sum), first(_first), last(_last), max_diff(_max_diff) {}

    bool empty() const {
        return max_diff < 0;
    }

    void apply(int length, const segment_change &change) {
        if (change.has_set()) {
            maximum = change.to_set;
            sum = int64_t(length) * change.to_set;
            first = last = change.to_set;
            max_diff = 0;
        }

        maximum += change.to_add;
        sum += int64_t(length) * change.to_add;
        first += change.to_add;
        last += change.to_add;
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
        max_diff = max({max_diff, other.max_diff, abs(last - other.first)});
        last = other.last;
    }

    // TODO: decide whether to re-implement this for better performance. Mainly relevant when segments contain arrays.
    void join(const segment &seg0, const segment &seg1) {
        *this = seg0;
        join(seg1);
    }
};

// TODO: persistent trees have both bad constant factor and high memory usage. Make sure we won't TLE or MLE.
struct persistent_seg_tree {
    // Note: when setting LAZY_PROPAGATION = false, use persistent_basic_seg_tree instead.
    static const bool LAZY_PROPAGATION = true;

    struct node {
        segment seg;
        segment_change change;
        int left = -1, right = -1;
    };

    int tree_n = 0;
    vector<node> tree;
    int reserve_size = INT32_MAX;

    persistent_seg_tree(int n = -1, int max_updates = 0) {
        if (n >= 0)
            init(n, max_updates);
    }

    void init(int n, int max_updates = 0) {
        tree_n = 1;

        while (tree_n < n)
            tree_n *= 2;

        reserve_size = INT32_MAX;

        if (max_updates > 0) {
            int tree_height = 32 - __builtin_clz(tree_n);
            reserve_size = 2 * tree_n + max_updates * 4 * tree_height;
            tree.reserve(reserve_size);
        }

        tree.assign(2 * tree_n, {});

        for (int i = 1; i < tree_n; i++) {
            tree[i].left = 2 * i;
            tree[i].right = 2 * i + 1;
        }
    }

    int left(int position) const { return tree[position].left; }
    int right(int position) const { return tree[position].right; }
    segment& seg(int position) { return tree[position].seg; }
    const segment& seg(int position) const { return tree[position].seg; }
    segment_change& change(int position) { return tree[position].change; }
    const segment_change& change(int position) const { return tree[position].change; }

    // Builds our tree from an array in O(n).
    void build(const vector<segment> &initial) {
        int n = int(initial.size());
        assert(n <= tree_n);

        for (int i = 0; i < n; i++)
            seg(tree_n + i) = initial[i];

        for (int position = tree_n - 1; position > 0; position--)
            seg(position).join(seg(left(position)), seg(right(position)));
    }

    int _make_copy(int position) {
        assert(0 <= position && position < int(tree.size()));
        tree.push_back(tree[position]);
        assert(int(tree.size()) <= reserve_size);
        return int(tree.size()) - 1;
    }

    void _push_down(int position, int length) {
        if (!change(position).has_change())
            return;

        seg(left(position)).apply(length / 2, change(position));
        seg(right(position)).apply(length / 2, change(position));
        change(left(position)) = change(left(position)).combine(change(position));
        change(right(position)) = change(right(position)).combine(change(position));
        change(position) = segment_change();
    }

    void _query_tree(int position, int start, int end, int a, int b, segment &answer, const segment_change &propagate) const {
        if (a <= start && end <= b) {
            segment current = seg(position);
            current.apply(end - start, propagate);
            answer.join(current);
            return;
        }

        if (left(position) < 0 || right(position) < 0)
            return;

        int mid = (start + end) / 2;
        segment_change next_propagate = change(position).combine(propagate);
        if (a < mid) _query_tree(left(position), start, mid, a, b, answer, next_propagate);
        if (b > mid) _query_tree(right(position), mid, end, a, b, answer, next_propagate);
    }

    segment query(int root, int a, int b) const {
        assert(root > 0 && 0 <= a && a <= b && b <= tree_n);
        segment answer;
        _query_tree(root, 0, tree_n, a, b, answer, segment_change());
        return answer;
    }

    // Directly assigning `tree[position].left = _make_copy(left(position))` results in segmentation faults, because the
    // address for `tree[position]` can be computed before calling `_make_copy`, which may reallocate `tree`.
    void _set_left(int position, int result) { tree[position].left = result; }
    void _set_right(int position, int result) { tree[position].right = result; }

    int _update_tree(int position, int start, int end, int a, int b, const segment_change &request, bool needs_copy) {
        if (needs_copy)
            position = _make_copy(position);

        if (a <= start && end <= b) {
            seg(position).apply(end - start, request);
            change(position) = change(position).combine(request);
            return position;
        }

        if (left(position) < 0 || right(position) < 0)
            return position;

        int mid = (start + end) / 2;

        if (LAZY_PROPAGATION) {
            _set_left(position, _make_copy(left(position)));
            _set_right(position, _make_copy(right(position)));
            _push_down(position, end - start);
            needs_copy = false;
        }

        if (a < mid) _set_left(position, _update_tree(left(position), start, mid, a, b, request, needs_copy));
        if (b > mid) _set_right(position, _update_tree(right(position), mid, end, a, b, request, needs_copy));
        seg(position).join(seg(left(position)), seg(right(position)));
        return position;
    }

    int update(int root, int a, int b, const segment_change &change) {
        assert(root > 0 && 0 <= a && a <= b && b <= tree_n);
        assert(LAZY_PROPAGATION || b - a <= 1);
        return _update_tree(root, 0, tree_n, a, b, change, true);
    }

    // Removes all updates starting from `root` or later from the tree.
    void undo_updates(int root) {
        assert(root >= 2 * tree_n && root <= int(tree.size()));
        tree.resize(root);
    }

    vector<segment> to_array(int root) {
        assert(root > 0);
        vector<int> level = {root};
        int length = tree_n;

        while (left(level.front()) >= 0) {
            vector<int> new_level;
            new_level.reserve(2 * level.size());

            for (int x : level) {
                _push_down(x, length);
                new_level.push_back(left(x));
                new_level.push_back(right(x));
            }

            swap(level, new_level);
            length /= 2;
        }

        vector<segment> segs(level.size());

        for (int i = 0; i < int(level.size()); i++)
            segs[i] = seg(level[i]);

        return segs;
    }

    // Finds the end of the last subarray starting at `first` satisfying `should_join` via binary search in O(log n).
    template<typename T_bool>
    int find_last_subarray(int root, T_bool &&should_join, int n, int first = 0) const {
        assert(root > 0 && 0 <= first && first <= n);
        segment current;

        // Check the degenerate case.
        if (!should_join(current, current))
            return first - 1;

        return y_combinator([&](auto search, int position, int start, int end, const segment_change &propagate) -> int {
            if (end <= first) {
                return end;
            } else if (first <= start && end <= n) {
                segment candidate = seg(position);
                candidate.apply(end - start, propagate);

                if (should_join(current, candidate)) {
                    current.join(candidate);
                    return end;
                }
            }

            if (end - start == 1)
                return start;

            int mid = (start + end) / 2;
            segment_change next_propagate = change(position).combine(propagate);
            int x = search(left(position), start, mid, next_propagate);
            return x < mid ? x : search(right(position), mid, end, next_propagate);
        })(root, 0, tree_n, segment_change());
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, Q;
    cin >> N >> Q;
    persistent_seg_tree tree(N, Q);
    tree.build(vector<segment>(N, segment(0, 0, 0, 0, 0)));
    int root = 1;

    for (int q = 0; q < Q; q++) {
        // This can handle the following operations (described below with inclusive 1-based indexing):
        // 1) "add a b x": add x to numbers[a] through numbers[b].
        // 2) "set a b x": set all of numbers[a] through numbers[b] to x.
        // 3) "max a b": compute the maximum of numbers[a] through numbers[b].
        // 4) "sum a b": compute the sum of numbers[a] through numbers[b].
        // 5) "diff a b": compute the maximum difference between any two consecutive numbers in numbers[a] through
        //      numbers[b]. If a = b, print 0.
        // 6) "fmax a x": given a, find the first index i with a <= i <= N such that numbers[i] >= x. If no such i
        //      exists, print -1.
        // 7) "fsum a x": given a, find the first index i with a - 1 <= i <= N such that numbers[a] + ... + numbers[i]
        //      is >= x (sum = 0 when i = a - 1). If no such i exists, print -1.
        string type;
        int a, b;
        int64_t x;
        cin >> type;

        if (type == "fmax") {
            cin >> a >> x;
            a--;

            int index = tree.find_last_subarray(root, [&](const segment &, const segment &add) -> bool {
                return add.maximum < x;
            }, N, a);

            assert(index <= N);
            cout << (index < N ? index + 1 : -1) << '\n';
            continue;
        } else if (type == "fsum") {
            cin >> a >> x;
            a--;

            int index = tree.find_last_subarray(root, [&](const segment &current, const segment &add) -> bool {
                return current.sum + add.sum < x;
            }, N, a);

            assert(index <= N);
            cout << (index < N ? index + 1 : -1) << '\n';
            continue;
        }

        cin >> a >> b;
        a--;
        assert(0 <= a && a < b && b <= N);

        if (type == "add" || type == "set")
            cin >> x;

        if (type == "add") {
            root = tree.update(root, a, b, segment_change(int(x)));
        } else if (type == "set") {
            root = tree.update(root, a, b, segment_change(0, int(x)));
        } else {
            segment seg = tree.query(root, a, b);

            if (type == "max")
                cout << seg.maximum << '\n';
            else if (type == "sum")
                cout << seg.sum << '\n';
            else if (type == "diff")
                cout << seg.max_diff << '\n';
            else
                assert(false);
        }
    }

    vector<segment> segs = tree.to_array(root);

    for (int i = 0; i < N; i++)
        cout << segs[i].sum << (i < N - 1 ? ' ' : '\n');
}
