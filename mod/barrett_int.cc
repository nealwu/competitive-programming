#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

struct barrett_reduction {
    unsigned mod;
    uint64_t div;

    barrett_reduction(unsigned m) : mod(m), div(-1LLU / m) {}

    unsigned operator()(uint64_t a) const {
#ifdef __SIZEOF_INT128__
        uint64_t q = uint64_t(__uint128_t(div) * a >> 64);
        uint64_t r = a - q * mod;
        return unsigned(r < mod ? r : r - mod);
#endif
        return unsigned(a % mod);
    }
};

template<const int &MOD, const barrett_reduction &barrett>
struct _b_int {
    int val;

    _b_int(int64_t v = 0) {
        if (v < 0) v = v % MOD + MOD;
        if (v >= MOD) v %= MOD;
        val = int(v);
    }

    _b_int(uint64_t v) {
        if (v >= uint64_t(MOD)) v %= MOD;
        val = int(v);
    }

    _b_int(int v) : _b_int(int64_t(v)) {}
    _b_int(unsigned v) : _b_int(uint64_t(v)) {}

    static int inv_mod(int a, int m = MOD) {
        // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Example
        int g = m, r = a, x = 0, y = 1;

        while (r != 0) {
            int q = g / r;
            g %= r; swap(g, r);
            x -= q * y; swap(x, y);
        }

        return x < 0 ? x + m : x;
    }

    explicit operator int() const { return val; }
    explicit operator unsigned() const { return val; }
    explicit operator int64_t() const { return val; }
    explicit operator uint64_t() const { return val; }
    explicit operator double() const { return val; }
    explicit operator long double() const { return val; }

    _b_int& operator+=(const _b_int &other) {
        val -= MOD - other.val;
        if (val < 0) val += MOD;
        return *this;
    }

    _b_int& operator-=(const _b_int &other) {
        val -= other.val;
        if (val < 0) val += MOD;
        return *this;
    }

    static unsigned fast_mod(uint64_t x) {
#if !defined(_WIN32) || defined(_WIN64)
        return barrett(x);
#endif
        // Optimized mod for Codeforces 32-bit machines.
        // x must be less than 2^32 * MOD for this to work, so that x / MOD fits in an unsigned 32-bit int.
        unsigned x_high = unsigned(x >> 32), x_low = unsigned(x);
        unsigned quot, rem;
        asm("divl %4\n"
            : "=a" (quot), "=d" (rem)
            : "d" (x_high), "a" (x_low), "r" (MOD));
        return rem;
    }

    _b_int& operator*=(const _b_int &other) {
        val = fast_mod(uint64_t(val) * other.val);
        return *this;
    }

    _b_int& operator/=(const _b_int &other) {
        return *this *= other.inv();
    }

    friend _b_int operator+(const _b_int &a, const _b_int &b) { return _b_int(a) += b; }
    friend _b_int operator-(const _b_int &a, const _b_int &b) { return _b_int(a) -= b; }
    friend _b_int operator*(const _b_int &a, const _b_int &b) { return _b_int(a) *= b; }
    friend _b_int operator/(const _b_int &a, const _b_int &b) { return _b_int(a) /= b; }

    _b_int& operator++() {
        val = val == MOD - 1 ? 0 : val + 1;
        return *this;
    }

    _b_int& operator--() {
        val = val == 0 ? MOD - 1 : val - 1;
        return *this;
    }

    _b_int operator++(int) { _b_int before = *this; ++*this; return before; }
    _b_int operator--(int) { _b_int before = *this; --*this; return before; }

    _b_int operator-() const {
        return val == 0 ? 0 : MOD - val;
    }

    friend bool operator==(const _b_int &a, const _b_int &b) { return a.val == b.val; }
    friend bool operator!=(const _b_int &a, const _b_int &b) { return a.val != b.val; }
    friend bool operator<(const _b_int &a, const _b_int &b) { return a.val < b.val; }
    friend bool operator>(const _b_int &a, const _b_int &b) { return a.val > b.val; }
    friend bool operator<=(const _b_int &a, const _b_int &b) { return a.val <= b.val; }
    friend bool operator>=(const _b_int &a, const _b_int &b) { return a.val >= b.val; }

    _b_int inv() const {
        return inv_mod(val);
    }

    _b_int pow(int64_t p) const {
        if (p < 0)
            return inv().pow(-p);

        _b_int a = *this, result = 1;

        while (p > 0) {
            if (p & 1)
                result *= a;

            p >>= 1;

            if (p > 0)
                a *= a;
        }

        return result;
    }

    friend ostream& operator<<(ostream &os, const _b_int &m) {
        return os << m.val;
    }

    friend istream& operator>>(istream &is, _b_int &m) {
        int64_t x;
        is >> x;
        m = x;
        return is;
    }
};

// TODO: double check all pre-initialized values.
// TODO: remember to re-initialize the `barrett` object after reading in `MOD`.
int MOD = int(1e9) + 7;
barrett_reduction barrett(MOD);
using barrett_int = _b_int<MOD, barrett>;


#include <chrono>
#include <random>

uint64_t random_address() { char *p = new char; delete p; return uint64_t(p); }
mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count() * (random_address() | 1));

bool is_prime(int64_t n) {
    if (n < 2)
        return false;

    for (int64_t p = 2; p * p <= n; p += p % 2 + 1)
        if (n % p == 0)
            return false;

    return true;
}

int main() {
    int mod;
    cin >> mod;
    barrett = barrett_reduction(MOD = mod);

    if (is_prime(MOD)) {
        for (int iter = 0; iter < 1000; iter++) {
            int64_t r = uniform_int_distribution<int64_t>(1, MOD - 1)(rng);
            int64_t inv1 = int64_t(barrett_int(r).inv());
            int64_t inv2 = int64_t(barrett_int(r).pow(MOD - 2));
            assert(inv1 == inv2);
            assert(r * inv1 % MOD == 1);
        }
    }

    barrett_int a, b;
    cin >> a >> b;
    a += b;
    cout << a << '\n';
    cin >> a >> b;
    a -= b;
    cout << a << '\n';
    cin >> a >> b;
    a *= b;
    cout << a << '\n';
    cin >> a >> b;
    if (__gcd(int(b), MOD) % MOD == 1) {
        a /= b;
        cout << a << '\n';
    } else {
        cout << "bad" << '\n';
    }
    cin >> a >> b;
    cout << a + b << '\n';
    cin >> a >> b;
    cout << a - b << '\n';
    cin >> a >> b;
    cout << a * b << '\n';
    cin >> a >> b;
    if (__gcd(int(b), MOD) % MOD == 1)
        cout << a / b << '\n';
    else
        cout << "bad" << '\n';
    cin >> a >> b;
    cout << -a << ' ' << -b << '\n';
}
