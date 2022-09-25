#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>
using namespace std;

uint64_t random_address() { char *p = new char; delete p; return uint64_t(p); }

uint64_t splitmix64(uint64_t x) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
}

const uint64_t FIXED_RANDOM = splitmix64(chrono::steady_clock::now().time_since_epoch().count() * (random_address() | 1));

template<typename T>
struct array_hash {
    using hash_t = uint64_t;

    int n;
    vector<T> arr;
    hash_t hash;

    array_hash(int _n = 0) {
        init(_n);
    }

    array_hash(const vector<T> &_arr) {
        init(_arr);
    }

    hash_t get_hash(int index) const {
        assert(0 <= index && index < n);
        // This ties arr[index] closely with index and ensures we call splitmix64 in between any two arithmetic operations.
        return splitmix64(arr[index] ^ splitmix64(index ^ FIXED_RANDOM));
    }

    void compute_hash() {
        hash = 0;

        for (int i = 0; i < n; i++)
            hash += get_hash(i);
    }

    void init(int _n) {
        n = _n;
        arr.assign(n, 0);
        compute_hash();
    }

    void init(const vector<T> &_arr) {
        arr = _arr;
        n = int(arr.size());
        compute_hash();
    }

    const T& operator[](int index) const {
        return arr[index];
    }

    void modify(int index, const T &value) {
        hash -= get_hash(index);
        arr[index] = value;
        hash += get_hash(index);
    }
};


#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;

int main() {
    int N, Q;
    cin >> N >> Q;
    vector<int64_t> A(N);

    for (auto &a : A)
        cin >> a;

    array_hash hasher(A);
    gp_hash_table<uint64_t, null_type> hashes;
    hashes.insert(hasher.hash);

    for (int q = 0; q < Q; q++) {
        int index;
        int64_t value;
        cin >> index >> value;
        index--;
        hasher.modify(index, value);
        assert(hasher[index] == value);
        hashes.insert(hasher.hash);
        cout << hashes.size() << '\n';
    }
}
