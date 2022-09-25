#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

template<typename T, bool maximum_mode = false>
struct RMQ {
    static int highest_bit(unsigned x) {
        return x == 0 ? -1 : 31 - __builtin_clz(x);
    }

    int n = 0;
    vector<T> values;
    vector<vector<int>> range_low;

    RMQ(const vector<T> &_values = {}) {
        if (!_values.empty())
            build(_values);
    }

    // Note: when `values[a] == values[b]`, returns b.
    int better_index(int a, int b) const {
        return (maximum_mode ? values[b] < values[a] : values[a] < values[b]) ? a : b;
    }

    void build(const vector<T> &_values) {
        values = _values;
        n = int(values.size());
        int levels = highest_bit(n) + 1;
        range_low.resize(levels);

        for (int k = 0; k < levels; k++)
            range_low[k].resize(n - (1 << k) + 1);

        for (int i = 0; i < n; i++)
            range_low[0][i] = i;

        for (int k = 1; k < levels; k++)
            for (int i = 0; i <= n - (1 << k); i++)
                range_low[k][i] = better_index(range_low[k - 1][i], range_low[k - 1][i + (1 << (k - 1))]);
    }

    // Note: breaks ties by choosing the largest index.
    int query_index(int a, int b) const {
        assert(0 <= a && a < b && b <= n);
        int level = highest_bit(b - a);
        return better_index(range_low[level][a], range_low[level][b - (1 << level)]);
    }

    T query_value(int a, int b) const {
        return values[query_index(a, b)];
    }
};

template<typename T, bool maximum_mode = false>
struct block_RMQ {
    static int highest_bit(unsigned x) {
        return x == 0 ? -1 : 31 - __builtin_clz(x);
    }

    // TODO: when adjusting BLOCK, make sure to adjust mask_t to match (e.g., BLOCK = 16 -> mask_t = uint16_t).
    static const int BLOCK = 8;
    using mask_t = uint8_t;
    using index_t = uint8_t;

    int n = 0;
    vector<T> values;
    vector<index_t> block_index;
    vector<mask_t> block_mask;
    RMQ<T, maximum_mode> rmq;

    block_RMQ(const vector<T> &_values = {}) {
        if (!_values.empty())
            build(_values);
    }

    static int get_block_start(int index) {
        return index - index % BLOCK;
    }

    // Note: when `values[a] == values[b]` and `maximum_mode` is false, returns `b`.
    int better_index(int a, int b) const {
        return (maximum_mode ? values[b] < values[a] : values[a] < values[b]) ? a : b;
    }

    void build(const vector<T> &_values) {
        values = _values;
        n = int(values.size());
        block_mask.assign(n + 1, 0);
        int mask = 0;

        for (int i = 0; i < n; i++) {
            if (i % BLOCK == 0)
                mask = 0;

            int block_start = get_block_start(i);

            while (mask != 0 && better_index(block_start + highest_bit(mask), i) == i)
                mask -= 1 << highest_bit(mask);

            mask |= 1 << (i - block_start);
            block_mask[i + 1] = mask_t(mask);
        }

        block_index.resize(n / BLOCK);
        vector<T> block_values(n / BLOCK);

        for (int i = 0; i < n / BLOCK; i++) {
            int block_start = BLOCK * i, best = block_start;

            for (int j = block_start + 1; j < block_start + BLOCK; j++)
                best = better_index(best, j);

            block_values[i] = values[best];
            block_index[i] = index_t(best - block_start);
        }

        rmq.build(block_values);
    }

    // Note: breaks ties by choosing the largest index.
    int query_index(int a, int b) const {
        assert(0 <= a && a < b && b <= n);
        int answer = a;
        int a_block_start = get_block_start(a);

        // Check if we are inside a single block.
        if (a_block_start == get_block_start(b - 1))
            return a + __builtin_ctz(block_mask[b] >> (a - a_block_start));

        if (a != a_block_start)
            answer = a + __builtin_ctz(block_mask[a_block_start + BLOCK] >> (a - a_block_start));

        int a_block = (a + BLOCK - 1) / BLOCK, b_block = b / BLOCK;

        if (a_block < b_block) {
            int rmq_index = rmq.query_index(a_block, b_block);
            answer = better_index(answer, BLOCK * rmq_index + block_index[rmq_index]);
        }

        int b_block_start = get_block_start(b);

        if (b != b_block_start)
            answer = better_index(answer, b_block_start + __builtin_ctz(block_mask[b]));

        return answer;
    }

    T query_value(int a, int b) const {
        return values[query_index(a, b)];
    }
};


#include <chrono>
#include <ctime>
#include <iomanip>
#include <random>

template<typename A, typename B> ostream& operator<<(ostream &os, const pair<A, B> &p) { return os << '(' << p.first << ", " << p.second << ')'; }
template<typename T_container, typename T = typename enable_if<!is_same<T_container, string>::value, typename T_container::value_type>::type> ostream& operator<<(ostream &os, const T_container &v) { os << '{'; string sep; for (const T &x : v) os << sep << x, sep = ", "; return os << '}'; }

uint64_t random_address() { char *p = new char; delete p; return uint64_t(p); }
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count() * (random_address() | 1));

// Uniformly distributed real number in [a, b).
double real_rng(double a = 0, double b = 1) {
    assert(a <= b);
    return uniform_real_distribution<double>(a, b)(rng);
}

// Uniformly distributed integer in [a, b].
int64_t unif_rng(int64_t a, int64_t b) {
    assert(a <= b);
    return uniform_int_distribution<int64_t>(a, b)(rng);
}

// Log-uniform distributed integer in [a, b]. P(a) > P(a + 1) > ... > P(b).
int64_t log_rng(int64_t a, int64_t b) {
    assert(a <= b);
    double min_val = double(a) - 0.5, max_val = double(b) + 0.5;
    int64_t x = int64_t(round(min_val - 1 + exp(real_rng(0, log(max_val - min_val + 1)))));

    // If x - a is large, randomize the lower bits in order to make up for double imprecision.
    static const int UNCHANGED_BITS = 30;

    if (uint64_t(x - a) >= 1LLU << UNCHANGED_BITS)
        x ^= rng() >> (__builtin_clzll(x - a) + UNCHANGED_BITS);

    return min(max(x, a), b);
}

// Returns +1 or -1 with 50% probability.
int sign_rng() {
    return 2 * int(unif_rng(0, 1)) - 1;
}


bool log_mode = false;

vector<string> process_args(int argc, char **argv) {
    vector<string> args;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];

        if (arg.find("-log") != string::npos)
            log_mode = true;
        else
            args.push_back(arg);
    }

    return args;
}

void do_queries() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N, Q;
    cin >> N >> Q;
    vector<int64_t> values(N);

    for (int64_t &v : values)
        cin >> v;

    block_RMQ<int64_t> rmq_min(values);
    block_RMQ<int64_t, true> rmq_max(values);

    for (int q = 0; q < Q; q++) {
        int a, b;
        cin >> a >> b;
        a--;
        cout << rmq_min.query_value(a, b) << ' ' << rmq_max.query_value(a, b) << '\n';
    }
}

int main(int argc, char **argv) {
    vector<string> args = process_args(argc, argv);

    if (args.empty()) {
        do_queries();
        return 0;
    }

    int N = args.size() >= 1 ? stoi(args[0]) : int(1e6 + 3);
    int Q = args.size() >= 2 ? stoi(args[1]) : int(2e6);

    vector<int64_t> values(N);
    vector<pair<int, int>> queries(Q);

    for (int i = 0; i < N; i++)
        values[i] = (i == 0 || sign_rng() > 0) ? rng() : values[unif_rng(0, i - 1)];

    const int BLOCK = block_RMQ<int64_t>::BLOCK;

    for (int i = 0; i < Q; i++) {
        int len = int(unif_rng(1, N));
        int start = int(unif_rng(0, N - len));

        if (BLOCK > 2 && N >= BLOCK && unif_rng(1, 100) <= 20) {
            // Create the pathological case for the block_RMQ algorithm where we have to iterate over the full range.
            len = (BLOCK - 2) - int(log_rng(0, BLOCK - 3));

            do {
                start = int(unif_rng(0, N - len));
            } while ((start + BLOCK - 1) / BLOCK != (start + len) / BLOCK + 1);
        }

        queries[i] = {start, start + len};
    }

    cout << fixed << setprecision(3);
    long double begin, build_time, query_time, combined;

    begin = clock();
    RMQ<int64_t> rmq_min(values);
    RMQ<int64_t, true> rmq_max(values);
    build_time = (clock() - begin) / CLOCKS_PER_SEC;

    begin = clock();
    uint64_t rmq_hash = 0;

    for (int i = 0; i < Q; i++) {
        rmq_hash = 37 * rmq_hash + rmq_min.query_index(queries[i].first, queries[i].second);
        rmq_hash = 37 * rmq_hash + rmq_max.query_index(queries[i].first, queries[i].second);
    }

    query_time = (clock() - begin) / CLOCKS_PER_SEC;
    combined = build_time + query_time;

    cout << "RMQ      : " << combined << "s (" << build_time << "s build, " << query_time << "s query), hash = " << rmq_hash << endl;

    begin = clock();
    block_RMQ<int64_t> block_rmq_min(values);
    block_RMQ<int64_t, true> block_rmq_max(values);
    build_time = (clock() - begin) / CLOCKS_PER_SEC;

    begin = clock();
    uint64_t block_rmq_hash = 0;

    for (int i = 0; i < Q; i++) {
        block_rmq_hash = 37 * block_rmq_hash + block_rmq_min.query_index(queries[i].first, queries[i].second);
        block_rmq_hash = 37 * block_rmq_hash + block_rmq_max.query_index(queries[i].first, queries[i].second);
    }

    query_time = (clock() - begin) / CLOCKS_PER_SEC;
    combined = build_time + query_time;

    cout << fixed << setprecision(3);
    cout << "block_RMQ: " << combined << "s (" << build_time << "s build, " << query_time << "s query), hash = " << block_rmq_hash << endl;
    assert(rmq_hash == block_rmq_hash);
}
