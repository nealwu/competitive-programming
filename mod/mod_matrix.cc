#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
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

        do {
            product *= MOD - MOD / v;
            v = MOD % v;
        } while (v >= SAVE_INV);

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

const int MOD = 998244353;
using mod_int = _m_int<MOD>;


// TODO: if using mod_column_vector, we can write the mod_matrix in the format matrix[x] = a row of coefficients used to
// build the x-th element of the mod_column_vector. So matrix[0][2] is the coefficient that element 2 contributes to the
// next element 0.
// The other option is to take a single-row 1 * n mod_matrix and multiply it by the n * n mod_matrix. Then matrix[0][2]
// is the coefficient that 0 contributes to the next element 2.
struct mod_column_vector {
    int rows;
    vector<mod_int> values;

    mod_column_vector(int _rows = 0) {
        init(_rows);
    }

    template<typename T>
    mod_column_vector(const vector<T> &v) {
        init(v);
    }

    void init(int _rows) {
        rows = _rows;
        values.assign(rows, 0);
    }

    template<typename T>
    void init(const vector<T> &v) {
        rows = int(v.size());
        values = vector<mod_int>(v.begin(), v.end());
    }

    mod_int& operator[](int index) { return values[index]; }
    const mod_int& operator[](int index) const { return values[index]; }
};

// Warning: very inefficient for many small matrices of fixed size. For that, use mod_matrix_fixed_size.cc instead.
struct mod_matrix {
    static const uint64_t U64_BOUND = numeric_limits<uint64_t>::max() - uint64_t(MOD) * MOD;

    static mod_matrix IDENTITY(int n) {
        mod_matrix identity(n);

        for (int i = 0; i < n; i++)
            identity[i][i] = 1;

        return identity;
    }

    int rows, cols;
    vector<vector<mod_int>> values;

    mod_matrix(int _rows = 0, int _cols = -1) {
        init(_rows, _cols);
    }

    template<typename T>
    mod_matrix(const vector<vector<T>> &v) {
        init(v);
    }

    void init(int _rows, int _cols = -1) {
        rows = _rows;
        cols = _cols < 0 ? rows : _cols;
        values.assign(rows, vector<mod_int>(cols, 0));
    }

    template<typename T>
    void init(const vector<vector<T>> &v) {
        rows = int(v.size());
        cols = v.empty() ? 0 : int(v[0].size());
        values.assign(rows, vector<mod_int>(cols, 0));

        for (int i = 0; i < rows; i++) {
            assert(int(v[i].size()) == cols);
            copy(v[i].begin(), v[i].end(), values[i].begin());
        }
    }

    vector<mod_int>& operator[](int index) { return values[index]; }
    const vector<mod_int>& operator[](int index) const { return values[index]; }

    bool is_square() const {
        return rows == cols;
    }

    mod_matrix operator*(const mod_matrix &other) const {
        assert(cols == other.rows);
        mod_matrix product(rows, other.cols);
        vector<uint64_t> row;

        for (int i = 0; i < rows; i++) {
            row.assign(other.cols, 0);

            for (int j = 0; j < cols; j++)
                if (values[i][j] != 0)
                    for (int k = 0; k < other.cols; k++) {
                        row[k] += uint64_t(values[i][j]) * uint64_t(other[j][k]);

                        if (row[k] > U64_BOUND)
                            row[k] %= MOD;
                    }

            for (int k = 0; k < other.cols; k++)
                product[i][k] = row[k];
        }

        return product;
    }

    mod_matrix& operator*=(const mod_matrix &other) {
        return *this = *this * other;
    }

    mod_column_vector operator*(const mod_column_vector &column) const {
        assert(cols == column.rows);
        mod_column_vector product(rows);

        for (int i = 0; i < rows; i++) {
            uint64_t result = 0;

            for (int j = 0; j < cols; j++) {
                result += uint64_t(values[i][j]) * uint64_t(column[j]);

                if (result > U64_BOUND)
                    result %= MOD;
            }

            product[i] = result;
        }

        return product;
    }

    mod_matrix& operator*=(mod_int mult) {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i][j] *= mult;

        return *this;
    }

    mod_matrix operator*(mod_int mult) const {
        return mod_matrix(*this) *= mult;
    }

    mod_matrix& operator+=(const mod_matrix &other) {
        assert(rows == other.rows && cols == other.cols);

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i][j] += other[i][j];

        return *this;
    }

    mod_matrix operator+(const mod_matrix &other) const {
        return mod_matrix(*this) += other;
    }

    mod_matrix& operator-=(const mod_matrix &other) {
        assert(rows == other.rows && cols == other.cols);

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                values[i][j] -= other[i][j];

        return *this;
    }

    mod_matrix operator-(const mod_matrix &other) const {
        return mod_matrix(*this) -= other;
    }

    mod_matrix pow(int64_t p) const {
        assert(p >= 0);
        assert(is_square());
        mod_matrix m = *this, result = IDENTITY(rows);

        while (p > 0) {
            if (p & 1)
                result *= m;

            p >>= 1;

            if (p > 0)
                m *= m;
        }

        return result;
    }

    void print(ostream &os) const {
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                os << values[i][j] << (j < cols - 1 ? ' ' : '\n');

        os << '\n';
    }
};


void read_matrix(mod_matrix &m) {
    int r, c;
    cin >> r >> c;
    m = mod_matrix(r, c);

    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) {
            int x;
            cin >> x;
            m[i][j] = x;
        }
}

int main() {
    mod_matrix m1, m2;
    read_matrix(m1);
    read_matrix(m2);
    (m1 + m1).print(cout);
    (m2 - m2).print(cout);
    (m1 * m2).print(cout);

    read_matrix(m1);
    int64_t p;
    cin >> p;
    (m1 * mod_int(p)).print(cout);
    m1.pow(p).print(cout);
}
