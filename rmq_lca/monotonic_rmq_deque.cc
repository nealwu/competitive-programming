#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <queue>
#include <vector>
using namespace std;

const int INF = int(1e9) + 5;

template<typename T, bool maximum_mode, T T_INF>
struct monotonic_RMQ {
    deque<pair<T, int>> values;
    int current_index = 0;

    static bool is_better(T a, T b) {
        return maximum_mode ? b < a : a < b;
    }

    int prev_add_index = -INF, prev_query_index = -INF;

    // Adds a value and returns its index.
    int add(const T &x, int index = -INF) {
        if (index == -INF)
            index = current_index++;

        assert(index >= prev_add_index);
        prev_add_index = index;

        while (!values.empty() && !is_better(values.back().first, x))
            values.pop_back();

        values.emplace_back(x, index);
        return index;
    }

    // Queries for the maximum (minimum) with index at least the given `index`. `index` must be non-decreasing.
    // TODO: when calling query, make sure to handle the empty case.
    T query_index(int index) {
        assert(index >= prev_query_index);
        prev_query_index = index;

        while (!values.empty() && values.front().second < index)
            values.pop_front();

        if (values.empty())
            return maximum_mode ? (is_signed<T>::value ? -T_INF : 0) : T_INF;

        return values.front().first;
    }

    // Warning: does not work when adding with custom indices.
    T query_count(int count) {
        return query_index(current_index - count);
    }

    // Returns whether or not the queue has a result for querying the given index.
    bool has_index(int index) const {
        return !values.empty() && values.back().second >= index;
    }

    bool has_count(int count) const {
        return has_index(current_index - count);
    }
};


template<bool maximum_mode, typename T>
vector<T> rmq_every_k(const vector<T> &values, int k) {
    int n = int(values.size());
    assert(1 <= k && k <= n);
    monotonic_RMQ<T, maximum_mode, numeric_limits<T>::max()> rmq;
    vector<T> result(n - k + 1);

    for (int i = 0; i < n; i++) {
        rmq.add(values[i]);
        assert(rmq.has_index(i) && !rmq.has_index(i + 1));

        if (i >= k - 1) {
            assert(rmq.has_count(k));
            result[i - (k - 1)] = rmq.query_count(k);
        }
    }

    return result;
}


int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, K;
    cin >> N >> K;
    vector<int64_t> values(N);

    for (int64_t &val : values)
        cin >> val;

    vector<int64_t> min_result = rmq_every_k<false>(values, K);
    vector<int64_t> max_result = rmq_every_k<true>(values, K);

    for (int i = 0; i < int(min_result.size()); i++)
        cout << min_result[i] << ' ' << max_result[i] << '\n';
}
