#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
using namespace std;

// Enables online insertion of (key, value) pairs and querying of maximum value out of keys less than a given limit.
// To query minimums instead, set maximum_mode = false.
template<typename T_key, typename T_value, bool maximum_mode, T_value V_INF>
struct online_prefix_max {
    static bool is_better(T_value a, T_value b) {
        return maximum_mode ? b < a : a < b;
    }

    static T_value default_value() {
        return maximum_mode ? (is_signed<T_value>::value ? -V_INF : 0) : V_INF;
    }

    map<T_key, T_value> optimal;

    int size() const {
        return int(optimal.size());
    }

    // Queries the maximum value in the map over all entries with key < `key_limit`.
    // TODO: when calling query, make sure to handle the empty case.
    T_value query(T_key key_limit) const {
        auto it = optimal.lower_bound(key_limit);
        return it == optimal.begin() ? default_value() : prev(it)->second;
    }

    // Adds an entry to the map and discards entries that are now obsolete.
    void insert(T_key key, T_value value) {
        auto it = optimal.upper_bound(key);

        // Quit if value is suboptimal.
        if (it != optimal.begin() && !is_better(value, prev(it)->second))
            return;

        if (it != optimal.begin() && prev(it)->first == key)
            optimal.erase(prev(it));

        while (it != optimal.end() && !is_better(it->second, value))
            it = optimal.erase(it);

        optimal.insert(it, {key, value});
    }
};


template<typename T_online_prefix_max>
void merge_into(T_online_prefix_max &x, T_online_prefix_max &y) {
    if (x.size() < y.size())
        swap(x, y);

    // Note: merge in linear time when x and y are close in size to improve worst-case runtime by a factor of log log n:
    // n log^2 n -> n log^2 n / log log n
    for (auto &p : y.optimal)
        x.insert(p.first, p.second);

    y.optimal.clear();
}


const int64_t INF64 = int64_t(2e18) + 5;

int main() {
    cerr << fixed << setprecision(4);

    int N;
    cin >> N;
    vector<pair<int64_t, int64_t>> inputs(N);

    for (auto &input : inputs)
        cin >> input.first >> input.second;

    vector<int64_t> outputs;
    outputs.reserve(N);

    long double begin = clock();
    online_prefix_max<int64_t, int64_t, true, INF64> prefix_max;
    int max_size = 0;

    for (auto &input : inputs) {
        int64_t key = input.first, value = input.second;
        outputs.push_back(prefix_max.query(key));
        prefix_max.insert(key, value);
        max_size = max(max_size, prefix_max.size());
    }

    cerr << "max size = " << max_size << endl;
    cerr << (clock() - begin) / CLOCKS_PER_SEC << 's' << endl;

    for (auto &output : outputs)
        cout << output << '\n';
}
