#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

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

// Returns a number that is a1 mod m1 and a2 mod m2. Assumes m1 and m2 are relatively prime.
int64_t chinese_remainder_theorem(int64_t a1, int64_t m1, int64_t a2, int64_t m2) {
    if (m1 < m2)
        return chinese_remainder_theorem(a2, m2, a1, m1);

    // assert(__gcd(m1, m2) == 1);
    assert(m1 >= m2);
    int64_t k = (a2 - a1) % m2 * inv_mod(m1, m2) % m2;
    int64_t result = a1 + k * m1;

    if (result < 0)
        result += m1 * m2;

    assert(0 <= result && result < m1 * m2);
    assert(result % m1 == a1 && result % m2 == a2);
    return result;
}

template<typename T>
int64_t chinese_remainder_theorem(const vector<T> &a, const vector<T> &m) {
    assert(a.size() == m.size());
    int64_t result = a.front();
    int64_t mod = m.front();

    for (int i = 1; i < int(m.size()); i++) {
        result = chinese_remainder_theorem(result, mod, a[i], m[i]);
        mod *= m[i];
    }

    return result;
}


#include <chrono>
#include <random>

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

int main() {
    const int64_t LIMIT = 9e18;

    for (int64_t test = 1; ; test++) {
        int64_t m1, m2;

        do {
            m1 = log_rng(1, LIMIT);
            m2 = log_rng(1, LIMIT);
        } while ((long double) m1 * m2 > LIMIT || __gcd(m1, m2) != 1);

        int64_t a1 = rng() % m1;
        int64_t a2 = rng() % m2;
        int64_t result = chinese_remainder_theorem(a1, m1, a2, m2);
        assert(0 <= result && result < m1 * m2);
        assert(result % m1 == a1 && result % m2 == a2);
        assert(result == chinese_remainder_theorem(vector<int64_t>{a1, a2}, vector<int64_t>{m1, m2}));

        if (test % 1000000 == 0)
            cerr << test << " tests passed!" << endl;
    }
}
