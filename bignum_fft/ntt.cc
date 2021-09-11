#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <queue>
#include <vector>
using namespace std;

template<const int &MOD>
struct _m_int {
    int val;

    _m_int(int64_t v = 0) {
        if (v < 0) v = v % MOD + MOD;
        if (v >= MOD) v %= MOD;
        val = int(v);
    }

    _m_int(uint64_t v) {
        if (v >= MOD) v %= MOD;
        val = int(v);
    }

    _m_int(int v) : _m_int(int64_t(v)) {}
    _m_int(unsigned v) : _m_int(uint64_t(v)) {}

    explicit operator int() const { return val; }
    explicit operator unsigned() const { return val; }
    explicit operator int64_t() const { return val; }
    explicit operator uint64_t() const { return val; }
    explicit operator double() const { return val; }
    explicit operator long double() const { return val; }

    _m_int& operator+=(const _m_int &other) {
        val -= MOD - other.val;
        if (val < 0) val += MOD;
        return *this;
    }

    _m_int& operator-=(const _m_int &other) {
        val -= other.val;
        if (val < 0) val += MOD;
        return *this;
    }

    static unsigned fast_mod(uint64_t x, unsigned m = MOD) {
#if !defined(_WIN32) || defined(_WIN64)
        return unsigned(x % m);
#endif
        // Optimized mod for Codeforces 32-bit machines.
        // x must be less than 2^32 * m for this to work, so that x / m fits in an unsigned 32-bit int.
        unsigned x_high = unsigned(x >> 32), x_low = unsigned(x);
        unsigned quot, rem;
        asm("divl %4\n"
            : "=a" (quot), "=d" (rem)
            : "d" (x_high), "a" (x_low), "r" (m));
        return rem;
    }

    _m_int& operator*=(const _m_int &other) {
        val = fast_mod(uint64_t(val) * other.val);
        return *this;
    }

    _m_int& operator/=(const _m_int &other) {
        return *this *= other.inv();
    }

    friend _m_int operator+(const _m_int &a, const _m_int &b) { return _m_int(a) += b; }
    friend _m_int operator-(const _m_int &a, const _m_int &b) { return _m_int(a) -= b; }
    friend _m_int operator*(const _m_int &a, const _m_int &b) { return _m_int(a) *= b; }
    friend _m_int operator/(const _m_int &a, const _m_int &b) { return _m_int(a) /= b; }

    _m_int& operator++() {
        val = val == MOD - 1 ? 0 : val + 1;
        return *this;
    }

    _m_int& operator--() {
        val = val == 0 ? MOD - 1 : val - 1;
        return *this;
    }

    _m_int operator++(int) { _m_int before = *this; ++*this; return before; }
    _m_int operator--(int) { _m_int before = *this; --*this; return before; }

    _m_int operator-() const {
        return val == 0 ? 0 : MOD - val;
    }

    friend bool operator==(const _m_int &a, const _m_int &b) { return a.val == b.val; }
    friend bool operator!=(const _m_int &a, const _m_int &b) { return a.val != b.val; }
    friend bool operator<(const _m_int &a, const _m_int &b) { return a.val < b.val; }
    friend bool operator>(const _m_int &a, const _m_int &b) { return a.val > b.val; }
    friend bool operator<=(const _m_int &a, const _m_int &b) { return a.val <= b.val; }
    friend bool operator>=(const _m_int &a, const _m_int &b) { return a.val >= b.val; }

    static const int SAVE_INV = int(1e6) + 5;
    static _m_int save_inv[SAVE_INV];

    static void prepare_inv() {
        // Ensures that MOD is prime, which is necessary for the inverse algorithm below.
        for (int64_t p = 2; p * p <= MOD; p += p % 2 + 1)
            assert(MOD % p != 0);

        save_inv[0] = 0;
        save_inv[1] = 1;

        for (int i = 2; i < SAVE_INV; i++)
            save_inv[i] = save_inv[MOD % i] * (MOD - MOD / i);
    }

    _m_int inv() const {
        if (save_inv[1] == 0)
            prepare_inv();

        if (val < SAVE_INV)
            return save_inv[val];

        _m_int product = 1;
        int v = val;

        while (v >= SAVE_INV) {
            product *= MOD - MOD / v;
            v = MOD % v;
        }

        return product * save_inv[v];
    }

    _m_int pow(int64_t p) const {
        if (p < 0)
            return inv().pow(-p);

        _m_int a = *this, result = 1;

        while (p > 0) {
            if (p & 1)
                result *= a;

            p >>= 1;

            if (p > 0)
                a *= a;
        }

        return result;
    }

    friend ostream& operator<<(ostream &os, const _m_int &m) {
        return os << m.val;
    }
};

template<const int &MOD> _m_int<MOD> _m_int<MOD>::save_inv[_m_int<MOD>::SAVE_INV];

extern const int MOD = 998244353;
using mod_int = _m_int<MOD>;


template<const int &MOD>
struct NTT {
    using ntt_int = _m_int<MOD>;

    static int highest_bit(unsigned x) {
        return x == 0 ? -1 : 31 - __builtin_clz(x);
    }

    static bool is_power_of_two(int n) {
        return (n & (n - 1)) == 0;
    }

    static int round_up_power_two(int n) {
        while (n & (n - 1))
            n = (n | (n - 1)) + 1;

        return max(n, 1);
    }

    // Given n (a power of two), finds k such that n == 1 << k.
    static int get_length(int n) {
        assert(is_power_of_two(n));
        return __builtin_ctz(n);
    }

    vector<ntt_int> roots = {0, 1};
    vector<int> bit_reverse;
    int max_size = -1;
    ntt_int root;

    void reset() {
        roots = {0, 1};
        max_size = -1;
    }

    // Rearranges the indices to be sorted by lowest bit first, then second lowest, etc., rather than highest bit first.
    // This makes even-odd div-conquer much easier.
    void bit_reorder(int n, vector<ntt_int> &values) {
        if (int(bit_reverse.size()) != n) {
            bit_reverse.assign(n, 0);
            int length = get_length(n);

            for (int i = 1; i < n; i++)
                bit_reverse[i] = (bit_reverse[i >> 1] >> 1) | ((i & 1) << (length - 1));
        }

        for (int i = 0; i < n; i++)
            if (i < bit_reverse[i])
                swap(values[i], values[bit_reverse[i]]);
    }

    void find_root() {
        max_size = 1 << __builtin_ctz(MOD - 1);
        root = 2;

        // Find a max_size-th primitive root of MOD.
        while (!(root.pow(max_size) == 1 && root.pow(max_size / 2) != 1))
            root++;
    }

    void prepare_roots(int n) {
        if (max_size < 0)
            find_root();

        assert(n <= max_size);

        if (int(roots.size()) >= n)
            return;

        int length = get_length(int(roots.size()));
        roots.resize(n);

        // The roots array is set up such that for a given power of two n >= 2, roots[n / 2] through roots[n - 1] are
        // the first half of the n-th primitive roots of MOD.
        while (1 << length < n) {
            // z is a 2^(length + 1)-th primitive root of MOD.
            ntt_int z = root.pow(max_size >> (length + 1));

            for (int i = 1 << (length - 1); i < 1 << length; i++) {
                roots[2 * i] = roots[i];
                roots[2 * i + 1] = roots[i] * z;
            }

            length++;
        }
    }

    void fft_iterative(int n, vector<ntt_int> &values) {
        assert(is_power_of_two(n));
        prepare_roots(n);
        bit_reorder(n, values);

        for (int len = 1; len < n; len *= 2)
            for (int start = 0; start < n; start += 2 * len)
                for (int i = 0; i < len; i++) {
                    ntt_int even = values[start + i];
                    ntt_int odd = values[start + len + i] * roots[len + i];
                    values[start + len + i] = even - odd;
                    values[start + i] = even + odd;
                }
    }

    void invert_fft(int n, vector<ntt_int> &values) {
        ntt_int inv_n = ntt_int(n).inv();

        for (int i = 0; i < n; i++)
            values[i] *= inv_n;

        reverse(values.begin() + 1, values.end());
        fft_iterative(n, values);
    }

    // Note: `circular = true` can be used for a 2x speedup when only the `max(n, m) - min(n, m) + 1` fully overlapping
    // ranges are needed. It computes results using indices modulo the power-of-two FFT size; see the brute force below.
    template<typename T>
    vector<T> mod_multiply(const vector<T> &_left, const vector<T> &_right, bool circular = false) {
        if (_left.empty() || _right.empty())
            return {};

        vector<ntt_int> left(_left.begin(), _left.end());
        vector<ntt_int> right(_right.begin(), _right.end());

        int n = int(left.size());
        int m = int(right.size());

        int output_size = circular ? round_up_power_two(max(n, m)) : n + m - 1;
        int N = round_up_power_two(output_size);

        double brute_force_cost = 1.25 * n * m;
        double ntt_cost = 3.0 * N * (get_length(N) + 3);

        if (brute_force_cost < ntt_cost) {
            auto &&mod_output_size = [&](int x) {
                return x < output_size ? x : x - output_size;
            };

            static const uint64_t U64_BOUND = numeric_limits<uint64_t>::max() - uint64_t(MOD) * MOD;
            vector<uint64_t> result(output_size, 0);

            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++) {
                    int index = mod_output_size(i + j);
                    result[index] += uint64_t(left[i]) * uint64_t(right[j]);

                    if (result[index] > U64_BOUND)
                        result[index] %= MOD;
                }

            for (uint64_t &x : result)
                if (x >= MOD)
                    x %= MOD;

            return vector<T>(result.begin(), result.end());
        }

        left.resize(N, 0);
        right.resize(N, 0);

        if (left == right) {
            fft_iterative(N, left);
            right = left;
        } else {
            fft_iterative(N, left);
            fft_iterative(N, right);
        }

        for (int i = 0; i < N; i++)
            left[i] *= right[i];

        invert_fft(N, left);
        return vector<T>(left.begin(), left.begin() + output_size);
    }

    template<typename T>
    vector<T> mod_power(const vector<T> &v, int exponent) {
        assert(exponent >= 0);
        vector<T> result = {1};

        if (exponent == 0)
            return result;

        for (int k = highest_bit(exponent); k >= 0; k--) {
            result = mod_multiply(result, result);

            if (exponent >> k & 1)
                result = mod_multiply(result, v);
        }

        return result;
    }

    // Multiplies many polynomials whose total degree is n in O(n log^2 n).
    template<typename T>
    vector<T> mod_multiply_all(const vector<vector<T>> &polynomials) {
        if (polynomials.empty())
            return {1};

        struct compare_size {
            bool operator()(const vector<T> &x, const vector<T> &y) {
                return x.size() > y.size();
            }
        };

        priority_queue<vector<T>, vector<vector<T>>, compare_size> pq(polynomials.begin(), polynomials.end());

        while (pq.size() > 1) {
            vector<T> a = pq.top(); pq.pop();
            vector<T> b = pq.top(); pq.pop();
            pq.push(mod_multiply(a, b));
        }

        return pq.top();
    }
};

NTT<MOD> ntt;



// Using the Chinese remainder theorem, does 3x NTT to perform `mod_multiply` for any mod.
namespace multi_ntt {
    int64_t inv_mod(int64_t a, int64_t m) {
        // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Example
        int64_t g = m, r = a, x = 0, y = 1;

        while (r != 0) {
            int64_t q = g / r;
            g %= r; swap(g, r);
            x -= q * y; swap(x, y);
        }

        assert(g == 1);
        assert(y == m || y == -m);
        return x < 0 ? x + m : x;
    }

    array<int64_t, 2> prepare_triple_crt(array<int, 3> m) {
        array<int64_t, 2> inv;
        inv[0] = inv_mod(m[0], m[1]);
        inv[1] = inv_mod(int64_t(m[0]) * m[1], m[2]);
        return inv;
    }

    // Returns the number in [0, m0 * m1) that is a0 mod m0 and a1 mod m1.
    int64_t chinese_remainder_theorem(int64_t a0, int64_t m0, int64_t a1, int64_t m1, int64_t inv_m0) {
        int64_t k = (a1 - a0) % m1 * inv_m0 % m1;
        int64_t result = a0 + k * m0;

        if (result < 0)
            result += m0 * m1;

        return result;
    }

    // Returns the number in [0, m[0] * m[1] * m[2]) that is a[0] mod m[0], a[1] mod m[1], and a[2] mod m[2], modulo `mod`.
    int triple_crt(array<int, 3> a, array<int, 3> m, array<int64_t, 2> inv, int mod) {
        int64_t m01 = int64_t(m[0]) * m[1];
        int64_t m01_mod = m01 % mod;
        int64_t a01 = chinese_remainder_theorem(a[0], m[0], a[1], m[1], inv[0]);

        int64_t k = (a[2] - a01) % m[2] * inv[1] % m[2];
        int64_t mod_result = (a01 + k * m01_mod) % mod;
        long double double_result = a01 + (long double) k * m01;

        if (double_result < 0)
            mod_result = (mod_result + m01_mod * m[2]) % mod;

        if (mod_result < 0)
            mod_result += mod;

        return int(mod_result);
    }

    extern const int MOD2 = 1711276033;
    extern const int MOD3 = 2113929217;
    const array<int64_t, 2> INV = prepare_triple_crt({MOD, MOD2, MOD3});
    const int64_t INV23 = inv_mod(MOD2, MOD3);

    NTT<MOD2> ntt2;
    NTT<MOD3> ntt3;

    // Mod multiply with any mod using three NTTs.
    template<typename T>
    vector<T> mod_multiply(const vector<T> &left, const vector<T> &right, int mod, bool circular = false) {
        vector<T> product1 = ntt.mod_multiply(left, right, circular);
        vector<T> product2 = ntt2.mod_multiply(left, right, circular);
        vector<T> product3 = ntt3.mod_multiply(left, right, circular);

        for (int i = 0; i < int(product1.size()); i++)
            product1[i] = triple_crt({product1[i], product2[i], product3[i]}, {MOD, MOD2, MOD3}, INV, mod);

        return product1;
    }

    // Multiply exactly (no mods) using two NTTs. Supports outputs up to 3.6e18.
    template<typename T_out, typename T_in>
    vector<T_out> multiply(const vector<T_in> &left, const vector<T_in> &right, bool circular = false) {
        vector<T_in> product2 = ntt2.mod_multiply(left, right, circular);
        vector<T_in> product3 = ntt3.mod_multiply(left, right, circular);
        vector<T_out> product(product2.size());

        for (int i = 0; i < int(product2.size()); i++)
            product[i] = chinese_remainder_theorem(product2[i], MOD2, product3[i], MOD3, INV23);

        return product;
    }
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    string task;
    cin >> task;
    assert(task == "mod_multiply");

    int n, m, mod;
    bool circular;
    cin >> n >> m >> mod >> circular;
    vector<int> left(n), right(m);

    for (int i = 0; i < n; i++)
        cin >> left[i];

    for (int i = 0; i < m; i++)
        cin >> right[i];

    const int TRIPLE_CUTOFF = 100000;
    vector<int> answer;

    if (mod == MOD) {
        answer = ntt.mod_multiply(left, right, circular);
    } else if (mod < TRIPLE_CUTOFF) {
        vector<int64_t> product = multi_ntt::multiply<int64_t>(left, right, circular);

        for (auto &x : product)
            answer.push_back(int(x % mod));
    } else {
        answer = multi_ntt::mod_multiply(left, right, mod, circular);
    }

    for (auto &x : answer)
        cout << x << '\n';
}
