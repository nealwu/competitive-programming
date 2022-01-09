// TODO: pragmas may help make this faster.
// #pragma GCC optimize("Ofast,unroll-loops,no-stack-protector,fast-math,inline")

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

// search_buckets provides two operations on an array:
// 1) set array[i] = x
// 2) count how many i in [start, end) satisfy array[i] < value
// Both operations take sqrt(n log n) time. Amazingly, because of the cache efficiency this is faster than the online
// (log n)^2 algorithm until n = 2-5 million.
template<typename T>
struct search_buckets {
    // values are just the values in order. buckets are sorted in segments of bucket_size (last segment may be smaller)
    int n, bucket_size;
    vector<T> values, buckets;

    search_buckets(const vector<T> &initial = {}) {
        init(initial);
    }

    int get_bucket_start(int index) const {
        return index - index % bucket_size;
    }

    int get_bucket_end_from_start(int bucket_start) const {
        return min(bucket_start + bucket_size, n);
    }

    void init(const vector<T> &initial) {
        values = buckets = initial;
        n = int(values.size());
        bucket_size = int(3 * sqrt(n * log(n + 1)) + 1);
        cerr << "Bucket size: " << bucket_size << endl;

        for (int start = 0; start < n; start += bucket_size)
            sort(buckets.begin() + start, buckets.begin() + get_bucket_end_from_start(start));
    }

    int bucket_count_less_than(int bucket_start, T value) const {
        auto begin = buckets.begin() + bucket_start;
        auto end = buckets.begin() + get_bucket_end_from_start(bucket_start);
        return int(lower_bound(begin, end, value) - begin);
    }

    int count_less_than(int start, int end, T value) const {
        int count = 0;
        int bucket_start = get_bucket_start(start);
        int bucket_end = get_bucket_end_from_start(bucket_start);

        if (start - bucket_start < bucket_end - start) {
            while (start > bucket_start)
                count -= values[--start] < value;
        } else {
            while (start < bucket_end)
                count += values[start++] < value;
        }

        bucket_start = get_bucket_start(end);
        bucket_end = get_bucket_end_from_start(bucket_start);

        if (end - bucket_start < bucket_end - end) {
            while (end > bucket_start)
                count += values[--end] < value;
        } else {
            while (end < bucket_end)
                count -= values[end++] < value;
        }

        while (start < end) {
            count += bucket_count_less_than(start, value);
            start = get_bucket_end_from_start(start);
        }

        assert(start == end);
        return count;
    }

    int prefix_count_less_than(int length, T value) const {
        return count_less_than(0, length, value);
    }

    void modify(int index, T value) {
        int bucket_start = get_bucket_start(index);
        int old_pos = bucket_start + bucket_count_less_than(bucket_start, values[index]);
        int new_pos = bucket_start + bucket_count_less_than(bucket_start, value);

        if (old_pos < new_pos) {
            copy(buckets.begin() + old_pos + 1, buckets.begin() + new_pos, buckets.begin() + old_pos);
            new_pos--;
            // memmove(&buckets[old_pos], &buckets[old_pos + 1], (new_pos - old_pos) * sizeof(T));
        } else {
            copy_backward(buckets.begin() + new_pos, buckets.begin() + old_pos, buckets.begin() + old_pos + 1);
            // memmove(&buckets[new_pos + 1], &buckets[new_pos], (old_pos - new_pos) * sizeof(T));
        }

        buckets[new_pos] = value;
        values[index] = value;
    }
};


// Solution to https://www.spoj.com/problems/RACETIME

int main() {
    int N, Q;
    scanf("%d %d", &N, &Q);
    vector<int> initial(N);

    for (int &a : initial)
        scanf("%d", &a);

    search_buckets<int> buckets(initial);

    for (int q = 0; q < Q; q++) {
        char type;
        int start, end, value;
        scanf(" %c", &type);

        if (type == 'M') {
            scanf("%d %d", &start, &value);
            start--;
            buckets.modify(start, value);
        } else if (type == 'C') {
            scanf("%d %d %d", &start, &end, &value);
            start--;
            printf("%d\n", buckets.count_less_than(start, end, value + 1));
        } else {
            assert(false);
        }
    }
}
