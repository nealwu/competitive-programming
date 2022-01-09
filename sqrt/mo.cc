// Solution to https://codeforces.com/contest/86/problem/D
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

using mo_value = int;
using mo_answer = int64_t;

// TODO: re-implement this struct.
struct mo_state {
    static const int A_MAX = int(1e6) + 5;

    int freq[A_MAX] = {0};
    int64_t sum = 0;

    void add_left(const mo_value &x) {
        sum += int64_t(2 * freq[x] + 1) * x;
        freq[x]++;
    }

    void add_right(const mo_value &x) {
        add_left(x);
    }

    void remove_left(const mo_value &x) {
        freq[x]--;
        sum -= int64_t(2 * freq[x] + 1) * x;
    }

    void remove_right(const mo_value &x) {
        remove_left(x);
    }

    mo_answer get_answer() const {
        return sum;
    }
};

struct mo_query {
    int start, end, block, index;

    mo_query() : start(0), end(0) {}

    mo_query(int _start, int _end) : start(_start), end(_end) {}

    bool operator<(const mo_query &other) const {
        if (block != other.block)
            return block < other.block;

        return block % 2 == 0 ? end < other.end : end > other.end;
    }
};

struct mo {
    int n, block_size;
    vector<mo_value> values;

    mo(vector<mo_value> initial = {}) {
        if (!initial.empty())
            init(initial);
    }

    void init(const vector<mo_value> &initial) {
        values = initial;
        n = int(values.size());
        block_size = int(1.6 * sqrt(n) + 1);
    }

    void update_state(mo_state &state, const mo_query &first, const mo_query &second) const {
        if (max(first.start, second.start) >= min(first.end, second.end)) {
            for (int i = first.start; i < first.end; i++)
                state.remove_left(values[i]);

            for (int i = second.start; i < second.end; i++)
                state.add_right(values[i]);

            return;
        }

        for (int i = first.start - 1; i >= second.start; i--)
            state.add_left(values[i]);

        for (int i = first.end; i < second.end; i++)
            state.add_right(values[i]);

        for (int i = first.start; i < second.start; i++)
            state.remove_left(values[i]);

        for (int i = first.end - 1; i >= second.end; i--)
            state.remove_right(values[i]);
    }

    vector<mo_answer> solve(vector<mo_query> queries) const {
        for (int i = 0; i < int(queries.size()); i++) {
            queries[i].index = i;
            queries[i].block = queries[i].start / block_size;
        }

        sort(queries.begin(), queries.end());
        mo_state state;
        mo_query last_query;
        vector<mo_answer> answers(queries.size());

        for (mo_query &q : queries) {
            update_state(state, last_query, q);
            answers[q.index] = state.get_answer();
            last_query = q;
        }

        return answers;
    }
};


int main() {
    int N, Q;
    scanf("%d %d", &N, &Q);
    vector<mo_value> A(N);

    for (mo_value &a : A)
        scanf("%d", &a);

    mo solver(A);
    vector<mo_query> queries(Q);

    for (mo_query &q : queries) {
        scanf("%d %d", &q.start, &q.end);
        q.start--;
    }

    vector<mo_answer> answers = solver.solve(queries);

    for (mo_answer &answer : answers)
        printf("%lld\n", (long long) answer);
}
