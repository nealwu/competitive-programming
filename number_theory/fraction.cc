#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;

struct fraction {
    // TODO: set this to false if it's unnecessary and the time limit might be tight.
    // CHECK_OVERFLOW64 = true can run up to 2 times slower (particularly on CF).
    static const bool CHECK_OVERFLOW64 = true;

    // TODO: consider setting AUTO_REDUCE = false for faster code. In this case, remember to call reduce() at the end.
    static const bool AUTO_REDUCE = true;

    static int cross_sign(const fraction &a, const fraction &b) {
        if (CHECK_OVERFLOW64) {
            long double double_value = (long double) a.numer * b.denom - (long double) b.numer * a.denom;

            if (abs(double_value) > 1e18)
                return double_value > 0 ? +1 : -1;
        }

        uint64_t uint64_value = (uint64_t) a.numer * b.denom - (uint64_t) b.numer * a.denom;
        int64_t actual = int64_t(uint64_value);
        return (actual > 0) - (actual < 0);
    }

    int64_t numer, denom;

    fraction(int64_t n = 0, int64_t d = 1) : numer(n), denom(d) {
        check_denom();

        if (AUTO_REDUCE)
            reduce();
    }

    void check_denom() {
        if (denom < 0) {
            numer = -numer;
            denom = -denom;
        }
    }

    void reduce() {
        int64_t g = __gcd(abs(numer), denom);
        numer /= g;
        denom /= g;
    }

    bool is_integer() const {
        return denom == 1 || (!AUTO_REDUCE && denom != 0 && numer % denom == 0);
    }

    friend fraction operator+(const fraction &a, const fraction &b) {
        return fraction(a.numer * b.denom + b.numer * a.denom, a.denom * b.denom);
    }

    friend fraction operator-(const fraction &a, const fraction &b) {
        return fraction(a.numer * b.denom - b.numer * a.denom, a.denom * b.denom);
    }

    friend fraction operator*(const fraction &a, const fraction &b) {
        return fraction(a.numer * b.numer, a.denom * b.denom);
    }

    friend fraction operator/(const fraction &a, const fraction &b) {
        return fraction(a.numer * b.denom, a.denom * b.numer);
    }

    fraction& operator+=(const fraction &other) { return *this = *this + other; }
    fraction& operator-=(const fraction &other) { return *this = *this - other; }
    fraction& operator*=(const fraction &other) { return *this = *this * other; }
    fraction& operator/=(const fraction &other) { return *this = *this / other; }

    fraction& operator++() { numer += denom; return *this; }
    fraction& operator--() { numer -= denom; return *this; }

    fraction operator++(int) { fraction before = *this; ++*this; return before; }
    fraction operator--(int) { fraction before = *this; --*this; return before; }

    fraction operator-() const {
        return fraction(-numer, denom);
    }

    fraction inv() const {
        return fraction(denom, numer);
    }

    bool operator==(const fraction &other) const { return cross_sign(*this, other) == 0; }
    bool operator!=(const fraction &other) const { return cross_sign(*this, other) != 0; }
    bool operator<(const fraction &other) const { return cross_sign(*this, other) < 0; }
    bool operator>(const fraction &other) const { return cross_sign(*this, other) > 0; }
    bool operator<=(const fraction &other) const { return cross_sign(*this, other) <= 0; }
    bool operator>=(const fraction &other) const { return cross_sign(*this, other) >= 0; }

    explicit operator double() const {
        return double(numer) / double(denom);
    }

    explicit operator long double() const {
        return (long double) numer / (long double) denom;
    }

    friend fraction abs(const fraction &f) {
        return fraction(abs(f.numer), f.denom);
    }

    friend ostream& operator<<(ostream& out, const fraction &frac) {
        return out << frac.numer << '/' << frac.denom;
    }
};


int main() {
    cout << setprecision(12);

    int64_t A, B, C, D;
    cin >> A >> B >> C >> D;
    fraction x(A, B), y(C, D);
    cout << min(x, y) << '\n';
    cout << max(x, y) << '\n';
    cout << x + y << '\n';
    cout << x - y << '\n';
    cout << x * y << '\n';
    cout << x / y << '\n';
    x++;
    y--;
    cout << x << '\n';
    cout << y << '\n';
    cout << abs(x) << '\n';
    cout << abs(y) << '\n';
    cout << double(x) << '\n';
    cout << (long double) y << '\n';
}
