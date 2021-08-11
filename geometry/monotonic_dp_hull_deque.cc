#include <algorithm>
#include <cassert>
#include <iostream>
#include <queue>
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

const int64_t INF64 = int64_t(2e18) + 5;

// monotonic_dp_hull enables you to do the following two operations in amortized O(1) time:
// 1. Insert a pair (a_i, b_i) into the structure. a_i must be non-decreasing.
// 2. For any value of x, query the maximum value of a_i * x + b_i. x must be non-decreasing.
// All values a_i, b_i, and x can be positive or negative.
struct monotonic_dp_hull {
    deque<point> points;

    void clear() {
        points.clear();
        prev_x = -INF64;
        prev_y = 1;
    }

    int size() const {
        return int(points.size());
    }

    void insert(const point &p) {
        assert(points.empty() || p.x >= points.back().x);

        if (!points.empty() && p.x == points.back().x) {
            if (p.y <= points.back().y)
                return;

            points.pop_back();
        }

        while (size() >= 2 && !left_turn(p, points.back(), points[size() - 2]))
            points.pop_back();

        points.push_back(p);
    }

    void insert(int64_t a, int64_t b) {
        insert(point(a, b));
    }

    int64_t prev_x = -INF64, prev_y = 1;

    // Queries the maximum value of ax + by.
    int64_t query(int64_t x, int64_t y = 1) {
        assert(size() > 0);
        assert(y > 0);
        assert(prev_x == -INF64 || x * prev_y >= prev_x * y);
        prev_x = x; prev_y = y;

        while (size() >= 2 && (points[1].x - points[0].x) * x + (points[1].y - points[0].y) * y >= 0)
            points.pop_front();

        return points[0].x * x + points[0].y * y;
    }
};


int main() {
    int Q;
    scanf("%d", &Q);
    monotonic_dp_hull hull;

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
}
