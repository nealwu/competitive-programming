#include <algorithm>
#include <cassert>
#include <iostream>
#include <set>
#include <vector>
using namespace std;

// TODO: set this to false if it's unnecessary and the time limit might be tight.
// CHECK_OVERFLOW64 = true can run up to 2 times slower (particularly on CF).
const bool CHECK_OVERFLOW64 = true;

struct point {
    int64_t x, y;

    point() : x(0), y(0) {}

    point(int64_t _x, int64_t _y) : x(_x), y(_y) {}

    point operator-(const point &other) const {
        return point(x - other.x, y - other.y);
    }
};

int cross_sign(const point &a, const point &b) {
    if (CHECK_OVERFLOW64) {
        long double double_value = (long double) a.x * b.y - (long double) b.x * a.y;

        if (abs(double_value) > 1e18)
            return double_value > 0 ? +1 : -1;
    }

    uint64_t uint64_value = (uint64_t) a.x * b.y - (uint64_t) b.x * a.y;
    int64_t actual = int64_t(uint64_value);
    return (actual > 0) - (actual < 0);
}

bool left_turn(const point &a, const point &b, const point &c) {
    return cross_sign(b - a, c - a) > 0;
}

// dp_hull enables you to do the following two operations in amortized O(log n) time:
// 1. Insert a pair (a_i, b_i) into the structure
// 2. For any value of x, query the maximum value of a_i * x + b_i
// All values a_i, b_i, and x can be positive or negative.
struct dp_hull {
    struct segment {
        point p;
        mutable point next_p;

        segment(point _p = {0, 0}) : p(_p), next_p(_p) {}

        // Note: this operator< supports `segments.lower_bound(point p)`. In order to support `upper_bound` as well, we
        // should also define `friend bool operator<(point p, segment s)`.
        bool operator<(const point &other) const {
            return (next_p.x - p.x) * other.x + (next_p.y - p.y) * other.y > 0;
        }

        bool operator<(const segment &other) const {
            return make_pair(p.x, p.y) < make_pair(other.p.x, other.p.y);
        }
    };

    set<segment, less<>> segments;

    int size() const {
        return int(segments.size());
    }

    bool bad(set<segment, less<>>::iterator it) const {
        if (it == segments.begin() || it == segments.end() || next(it) == segments.end())
            return false;

        return !left_turn(next(it)->p, it->p, prev(it)->p);
    }

    void insert(const point &p) {
        auto next_it = segments.lower_bound(segment(p));

        if (next_it != segments.end() && p.x == next_it->p.x)
            return;

        if (next_it != segments.begin()) {
            auto prev_it = prev(next_it);

            if (p.x == prev_it->p.x)
                segments.erase(prev_it);
            else if (next_it != segments.end() && !left_turn(next_it->p, p, prev_it->p))
                return;
        }

        auto it = segments.insert(next_it, segment(p));

        while (it != segments.begin() && bad(prev(it)))
            segments.erase(prev(it));

        while (bad(next(it)))
            segments.erase(next(it));

        if (it != segments.begin())
            prev(it)->next_p = it->p;

        it->next_p = next(it) != segments.end() ? next(it)->p : it->p;
    }

    void insert(int64_t a, int64_t b) {
        insert(point(a, b));
    }

    // Queries the maximum value of ax + by.
    int64_t query(int64_t x, int64_t y = 1) const {
        assert(size() > 0);
        assert(y > 0);
        auto it = segments.lower_bound(point(x, y));
        return it->p.x * x + it->p.y * y;
    }
};


int main() {
    int Q;
    scanf("%d", &Q);
    dp_hull hull;

    for (int q = 0; q < Q; q++) {
        int type;
        scanf("%d", &type);

        if (type == 1) {
            int a, b;
            scanf("%d %d", &a, &b);
            hull.insert(a, b);
        } else if (type == 2) {
            int x;
            scanf("%d", &x);
            printf("%lld\n", (long long) hull.query(2 * x, 2) / 2);
        } else {
            assert(false);
        }
    }

    cerr << "size: " << hull.size() << endl;
}
