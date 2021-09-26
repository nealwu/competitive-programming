#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <iostream>
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
    void join(const segment &a, const segment &b) {
        *this = a;
        join(b);
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

    template<typename T_range_op>
    void process_range(int position, int start, int end, int a, int b, bool needs_join, T_range_op &&range_op) {
        if (a <= start && end <= b) {
            range_op(position, end - start);
            return;
        }

        if (position >= tree_n)
            return;

        push_down(position, end - start);
        int mid = (start + end) / 2;
        if (a < mid) process_range(2 * position, start, mid, a, b, needs_join, range_op);
        if (b > mid) process_range(2 * position + 1, mid, end, a, b, needs_join, range_op);
        if (needs_join) tree[position].join(tree[2 * position], tree[2 * position + 1]);
    }

    segment query(int a, int b) {
        assert(0 <= a && a <= b && b <= tree_n);
        segment answer;

        process_range(1, 0, tree_n, a, b, false, [&](int position, int) {
            answer.join(tree[position]);
        });

        return answer;
    }

    segment query_full() const {
        return tree[1];
    }

    void update(int a, int b, const segment_change &change) {
        assert(0 <= a && a <= b && b <= tree_n);

        process_range(1, 0, tree_n, a, b, true, [&](int position, int length) {
            apply_and_combine(position, length, change);
        });
    }

    vector<segment> to_array() {
        for (int i = 1; i < tree_n; i++)
            push_down(i, tree_n >> highest_bit(i));

        vector<segment> segs(tree_n);

        for (int i = 0; i < tree_n; i++)
            segs[i] = tree[tree_n + i];

        return segs;
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

    // Finds the last subarray starting at `first` that satisifes `should_join` via binary search in O(log n).
    template<typename T_bool>
    int find_last_subarray(T_bool &&should_join, int n, int first = 0) {
        assert(0 <= first && first <= n);
        segment current;

        // Check the degenerate case.
        if (!should_join(current, current))
            return first - 1;

        return y_combinator([&](auto search, int position, int start, int end) -> int {
            if (end <= first) {
                return end;
            } else if (first <= start && end <= n && should_join(current, tree[position])) {
                current.join(tree[position]);
                return end;
            } else if (end - start == 1) {
                return start;
            }

            push_down(position, end - start);
            int mid = (start + end) / 2;
            int left = search(2 * position, start, mid);
            return left < mid ? left : search(2 * position + 1, mid, end);
        })(1, 0, tree_n);
    }
};


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, Q;
    cin >> N >> Q;
    seg_tree tree(N);
    tree.build(vector<segment>(N, segment(0, 0, 0, 0, 0)));

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

            int index = tree.find_last_subarray([&](const segment &, const segment &add) {
                return add.maximum < x;
            }, N, a);

            cout << (index < N ? index + 1 : -1) << '\n';
            continue;
        } else if (type == "fsum") {
            cin >> a >> x;
            a--;

            int index = tree.find_last_subarray([&](const segment &current, const segment &add) {
                return current.sum + add.sum < x;
            }, N, a);

            cout << (index < N ? index + 1 : -1) << '\n';
            continue;
        }

        cin >> a >> b;
        a--;
        assert(0 <= a && a < b && b <= N);

        if (type == "add" || type == "set")
            cin >> x;

        if (type == "add") {
            tree.update(a, b, segment_change(int(x)));
        } else if (type == "set") {
            if (b - a == 1)
                tree.update_single(a, segment(int(x), x, int(x), int(x), 0));
            else
                tree.update(a, b, segment_change(0, int(x)));
        } else {
            segment seg = a == 0 && b == N ? tree.query_full() : tree.query(a, b);

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

    vector<segment> segs = tree.to_array();

    for (int i = 0; i < N; i++)
        cout << segs[i].sum << (i < N - 1 ? ' ' : '\n');
}
