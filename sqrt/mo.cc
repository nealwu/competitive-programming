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
    const vector<mo_value> &values;

    vector<int> freq;
    int64_t sum;

    mo_state(const vector<mo_value> &_values) : values(_values) {
        int maximum = values.empty() ? 0 : *max_element(values.begin(), values.end());
        freq.assign(maximum + 1, 0);
        sum = 0;
    }

    void add_left(int index) {
        mo_value x = values[index];
        sum += int64_t(2 * freq[x] + 1) * x;
        freq[x]++;
    }

    void add_right(int index) {
        add_left(index);
    }

    void remove_left(int index) {
        mo_value x = values[index];
        freq[x]--;
        sum -= int64_t(2 * freq[x] + 1) * x;
    }

    void remove_right(int index) {
        remove_left(index);
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
    int n;
    vector<mo_value> values;

    mo(const vector<mo_value> &_values = {}) : values(_values) {
        n = int(values.size());
    }

    void update_state(mo_state &state, const mo_query &first, const mo_query &second) const {
        if (max(first.start, second.start) >= min(first.end, second.end)) {
            for (int i = first.start; i < first.end; i++)
                state.remove_left(i);

            for (int i = second.start; i < second.end; i++)
                state.add_right(i);

            return;
        }

        for (int i = first.start - 1; i >= second.start; i--)
            state.add_left(i);

        for (int i = first.end; i < second.end; i++)
            state.add_right(i);

        for (int i = first.start; i < second.start; i++)
            state.remove_left(i);

        for (int i = first.end - 1; i >= second.end; i--)
            state.remove_right(i);
    }

    vector<mo_answer> solve(vector<mo_query> queries) const {
        int block_size = int(1.5 * n / sqrt(max(int(queries.size()), 1)) + 1);

        for (int i = 0; i < int(queries.size()); i++) {
            queries[i].index = i;
            queries[i].block = queries[i].start / block_size;
        }

        sort(queries.begin(), queries.end());
        mo_state state(values);
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

    for (mo_query &qry : queries) {
        scanf("%d %d", &qry.start, &qry.end);
        qry.start--;
    }

    vector<mo_answer> answers = solver.solve(queries);

    for (mo_answer &answer : answers)
        printf("%lld\n", (long long) answer);
}
