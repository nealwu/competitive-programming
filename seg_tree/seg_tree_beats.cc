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


// Note: value_t must be able to handle sums of elements / sums of updates, not just individual elements.
using value_t = int64_t;

struct segment_change {
    // Use a sentinel value rather than a boolean to save significant memory (4-8 bytes per object).
    // SENTINEL is also conveniently the minimum possible value, which is a no-op for to_max.
    static const value_t SENTINEL = numeric_limits<value_t>::lowest();

    // Note that to_set goes first, then to_max, and then to_add.
    // TODO: check if these values can overflow int.
    value_t to_set, to_max, to_add;

    // TODO: make sure the default constructor is the identity segment_change.
    segment_change(value_t _to_add = 0, value_t _to_set = SENTINEL, value_t _to_max = SENTINEL)
        : to_set(_to_set), to_max(_to_max), to_add(_to_add) {}

    bool has_set() const { return to_set != SENTINEL; }
    bool has_max() const { return to_max != SENTINEL; }

    bool has_change() const {
        return has_set() || has_max() || to_add != 0;
    }

    // Return the combined result of applying this segment_change followed by `other`.
    // TODO: make sure to check for sentinel values.
    segment_change combine(const segment_change &other) const {
        if (other.has_set())
            return other;

        segment_change combined = *this;

        // to_set, to_max, to_add, other.to_max, other.to_add ->
        // to_set, max(to_max, other.to_max - to_add), to_add + other.to_add
        if (other.has_max())
            combined.to_max = max(to_max, other.to_max - to_add);

        combined.to_add += other.to_add;
        return combined;
    }
};

struct segment {
    // TODO: check if these values can overflow int.
    value_t minimum, second_min;
    int min_count, second_count;

    value_t maximum, sum;

    // TODO: make sure the default constructor is the identity segment.
    segment(value_t value = 0, int count = 0) {
        minimum = count != 0 ? value : numeric_limits<value_t>::max();
        min_count = count;
        second_min = numeric_limits<value_t>::max();
        second_count = 0;

        maximum = count != 0 ? value : numeric_limits<value_t>::lowest();
        sum = value_t(value) * count;
    }

    bool empty() const {
        return min_count == 0;
    }

    // apply returns true if we can stop here. If false, we need to continue updating down.
    bool apply(int length, const segment_change &change) {
        if (change.has_set()) {
            minimum = change.to_set;
            min_count = length;
            second_min = numeric_limits<value_t>::max();
            second_count = 0;

            maximum = change.to_set;
            sum = value_t(change.to_set) * length;
        }

        if (change.has_max() && change.to_max > minimum) {
            if (second_count != 0 && change.to_max >= second_min)
                return false;

            sum += value_t(change.to_max - minimum) * min_count;
            minimum = change.to_max;
            maximum = max(maximum, change.to_max);
        }

        if (change.to_add != 0) {
            minimum += change.to_add;
            if (second_count != 0) second_min += change.to_add;
            maximum += change.to_add;
            sum += value_t(change.to_add) * length;
        }

        return true;
    }

    void join(const segment &other) {
        if (empty()) {
            *this = other;
            return;
        } else if (other.empty()) {
            return;
        }

        auto update_second_min = [&](value_t cand_min, int cand_count) -> void {
            if (cand_min < second_min) {
                second_min = cand_min;
                second_count = cand_count;
            } else if (cand_min == second_min) {
                second_count += cand_count;
            }
        };

        if (minimum == other.minimum) {
            min_count += other.min_count;
            update_second_min(other.second_min, other.second_count);
        } else if (minimum < other.minimum) {
            update_second_min(other.minimum, other.min_count);
        } else {
            second_min = minimum;
            second_count = min_count;
            update_second_min(other.second_min, other.second_count);
            minimum = other.minimum;
            min_count = other.min_count;
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

struct seg_tree_beats {
    static int highest_bit(unsigned x) {
        return x == 0 ? -1 : 31 - __builtin_clz(x);
    }

    int tree_n = 0;
    vector<segment> tree;
    vector<segment_change> changes;

    seg_tree_beats(int n = -1) {
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

    bool apply_and_combine(int position, int length, const segment_change &change) {
        if (tree[position].apply(length, change)) {
            if (position < tree_n)
                changes[position] = changes[position].combine(change);

            return true;
        }

        return false;
    }

    void push_down(int position, int length) {
        if (changes[position].has_change()) {
            bool success = true;
            success &= apply_and_combine(2 * position, length / 2, changes[position]);
            success &= apply_and_combine(2 * position + 1, length / 2, changes[position]);
            assert(success);
            changes[position] = segment_change();
        }
    }

    template<typename T_range_op>
    void process_range(int position, int start, int end, int a, int b, bool needs_join, T_range_op &&range_op) {
        // range_op returns true if we can finish here; if false, we need to continue recursing down.
        if (a <= start && end <= b && range_op(position, end - start))
            return;

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

        process_range(1, 0, tree_n, a, b, false, [&](int position, int) -> bool {
            answer.join(tree[position]);
            return true;
        });

        return answer;
    }

    segment query_full() const {
        return tree[1];
    }

    void update(int a, int b, const segment_change &change) {
        assert(0 <= a && a <= b && b <= tree_n);

        process_range(1, 0, tree_n, a, b, true, [&](int position, int length) -> bool {
            return apply_and_combine(position, length, change);
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

    vector<segment> to_array() {
        for (int i = 1; i < tree_n; i++)
            push_down(i, tree_n >> highest_bit(i));

        vector<segment> segs(tree_n);

        for (int i = 0; i < tree_n; i++)
            segs[i] = tree[tree_n + i];

        return segs;
    }

    // Finds the end of the last subarray starting at `first` satisfying `should_join` via binary search in O(log n).
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
    seg_tree_beats tree(N);
    tree.build(vector<segment>(N, segment(0, 1)));

    for (int q = 0; q < Q; q++) {
        // This can handle the following operations (described below with inclusive 1-based indexing):
        // 1) "set_max a b x": all of numbers[a] through numbers[b] do num = max(num, x).
        // 2) "set a b x": set all of numbers[a] through numbers[b] to x.
        // 3) "sum a b": compute the sum of numbers[a] through numbers[b].
        // 4) "min a b": compute the min of numbers[a] through numbers[b].
        // 5) "max a b": compute the max of numbers[a] through numbers[b].
        // 6) "add a b x": add x to all of numbers[a] through numbers[b]. Including this operation changes the
        //    complexity from amortized O(log n) per operation to amortized O(log^2 n) per operation. However it still
        //    seems fast in practice.
        // 7) "fmax a x": given a, find the first index i with a <= i <= N such that numbers[i] >= x. If no such i
        //      exists, print -1.
        // 8) "fsum a x": given a, find the first index i with a - 1 <= i <= N such that numbers[a] + ... + numbers[i]
        //      is >= x (sum = 0 when i = a - 1). If no such i exists, print -1.
        string type;
        int a, b;
        value_t x;
        cin >> type;

        if (type == "fmax") {
            cin >> a >> x;
            a--;

            int index = tree.find_last_subarray([&](const segment &, const segment &add) -> bool {
                return add.maximum < x;
            }, N, a);

            cout << (index < N ? index + 1 : -1) << '\n';
            continue;
        } else if (type == "fsum") {
            cin >> a >> x;
            a--;

            int index = tree.find_last_subarray([&](const segment &current, const segment &add) -> bool {
                return current.sum + add.sum < x;
            }, N, a);

            cout << (index < N ? index + 1 : -1) << '\n';
            continue;
        }

        cin >> a >> b;
        a--;

        if (type == "set_max") {
            cin >> x;
            tree.update(a, b, segment_change(0, segment_change::SENTINEL, x));
        } else if (type == "set") {
            cin >> x;

            if (b - a == 1)
                tree.update_single(a, segment(x, 1));
            else
                tree.update(a, b, segment_change(0, x));
        } else if (type == "add") {
            cin >> x;
            tree.update(a, b, segment_change(x));
        } else {
            segment seg = a == 0 && b == N ? tree.query_full() : tree.query(a, b);

            if (type == "sum")
                cout << seg.sum << '\n';
            else if (type == "min")
                cout << seg.minimum << '\n';
            else if (type == "max")
                cout << seg.maximum << '\n';
            else
                assert(false);
        }
    }

    vector<segment> segs = tree.to_array();

    for (int i = 0; i < N; i++)
        cout << segs[i].sum << (i < N - 1 ? ' ' : '\n');
}
