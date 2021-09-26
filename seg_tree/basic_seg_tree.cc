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


// TODO: segment_change can be eliminated entirely in favor of just updating with a new segment instead.
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

    void apply(const segment_change &change) {
        if (change.has_set()) {
            maximum = first = last = change.to_set;
            sum = change.to_set;
            max_diff = 0;
        }

        maximum += change.to_add;
        sum += change.to_add;
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

int right_half[32];

struct basic_seg_tree {
    // TODO: POWER_OF_TWO_MODE is necessary in order to call query_full() or to binary search the tree.
    static const bool POWER_OF_TWO_MODE = true;

    int tree_n = 0;
    vector<segment> tree;

    basic_seg_tree(int n = -1) {
        if (n >= 0)
            init(n);
    }

    void init(int n) {
        if (POWER_OF_TWO_MODE) {
            tree_n = 1;

            while (tree_n < n)
                tree_n *= 2;
        } else {
            tree_n = n;
        }

        tree.assign(2 * tree_n, segment());
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

    segment query(int a, int b) const {
        assert(0 <= a && a <= b && b <= tree_n);
        segment answer;
        int r_size = 0;

        for (a += tree_n, b += tree_n; a < b; a /= 2, b /= 2) {
            if (a & 1)
                answer.join(tree[a++]);

            if (b & 1)
                right_half[r_size++] = --b;
        }

        for (int i = r_size - 1; i >= 0; i--)
            answer.join(tree[right_half[i]]);

        return answer;
    }

    segment query_full() const {
        assert(POWER_OF_TWO_MODE);
        return tree[1];
    }

    segment query_single(int index) const {
        assert(0 <= index && index < tree_n);
        return tree[tree_n + index];
    }

    void join_up(int position) {
        while (position > 1) {
            position /= 2;
            tree[position].join(tree[2 * position], tree[2 * position + 1]);
        }
    }

    void update(int index, const segment_change &change) {
        assert(0 <= index && index < tree_n);
        int position = tree_n + index;
        tree[position].apply(change);
        join_up(position);
    }

    void update(int index, const segment &seg) {
        assert(0 <= index && index < tree_n);
        int position = tree_n + index;
        tree[position] = seg;
        join_up(position);
    }

    // Finds the last subarray starting at `first` that satisifes `should_join` via binary search in O(log n).
    template<typename T_bool>
    int find_last_subarray(T_bool &&should_join, int n, int first = 0) const {
        assert(POWER_OF_TWO_MODE);
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
    basic_seg_tree tree(N);
    tree.build(vector<segment>(N, segment(0, 0, 0, 0, 0)));

    for (int q = 0; q < Q; q++) {
        // This can handle the following operations (described below with inclusive 1-based indexing):
        // 1) "add a b x": add x to numbers[a] through numbers[b]. Must have a = b.
        // 2) "set a b x": set all of numbers[a] through numbers[b] to x. Must have a = b.
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

            assert(index <= N);
            cout << (index < N ? index + 1 : -1) << '\n';
            continue;
        } else if (type == "fsum") {
            cin >> a >> x;
            a--;

            int index = tree.find_last_subarray([&](const segment &current, const segment &add) {
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
            assert(b - a == 1);
            tree.update(a, segment_change(int(x)));
        } else if (type == "set") {
            assert(b - a == 1);
            tree.update(a, segment_change(0, int(x)));
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

    for (int i = 0; i < N; i++)
        cout << tree.tree[tree.tree_n + i].sum << (i < N - 1 ? ' ' : '\n');
}
