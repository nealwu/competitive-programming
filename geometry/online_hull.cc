#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <vector>
using namespace std;

using coord_t = int;

// Maintains the upper hull for coordinates that can span a range of roughly [-1e9, 1e9].
struct upper_hull {
    map<coord_t, coord_t> points;
    int64_t area = 0;

    int size() const {
        return int(points.size());
    }

    int64_t cross(coord_t x1, coord_t y1, coord_t x2, coord_t y2) const {
        return (int64_t) x1 * y2 - (int64_t) x2 * y1;
    }

    // Note: all areas are signed and doubled. For triangles, clockwise (not counter-clockwise) is positive.

    int64_t trapezoid_area(coord_t x1, coord_t y1, coord_t x2, coord_t y2) const {
        return (int64_t) (x2 - x1) * (y1 + y2);
    }

    int64_t trapezoid_area(pair<coord_t, coord_t> p1, pair<coord_t, coord_t> p2) const {
        return trapezoid_area(p1.first, p1.second, p2.first, p2.second);
    }

    int64_t triangle_area(coord_t x1, coord_t y1, coord_t x2, coord_t y2, coord_t x3, coord_t y3) const {
        return cross(x2 - x1, y2 - y1, x2 - x3, y2 - y3);
    }

    int64_t triangle_area(pair<coord_t, coord_t> p1, pair<coord_t, coord_t> p2, pair<coord_t, coord_t> p3) const {
        return triangle_area(p1.first, p1.second, p2.first, p2.second, p3.first, p3.second);
    }

    // Gets the area that a point is responsible for; that is, how much area would be lost if the point were removed.
    int64_t point_area(map<coord_t, coord_t>::iterator it) {
        bool has_prev = it != points.begin();
        bool has_next = next(it) != points.end();

        if (has_prev && has_next)
            return triangle_area(*prev(it), *it, *next(it));

        int64_t sum = 0;

        if (has_prev)
            sum += trapezoid_area(*prev(it), *it);

        if (has_next)
            sum += trapezoid_area(*it, *next(it));

        return sum;
    }

    bool bad(map<coord_t, coord_t>::iterator it) {
        if (it == points.begin() || it == points.end() || next(it) == points.end())
            return false;

        // True if the three points form a left turn or line up straight.
        return point_area(it) <= 0;
    }

    void erase(map<coord_t, coord_t>::iterator it) {
        area -= point_area(it);
        points.erase(it);
    }

    bool insert(coord_t x, coord_t y) {
        if (points.find(x) != points.end()) {
            if (y <= points[x])
                return false;

            erase(points.find(x));
        }

        map<coord_t, coord_t>::iterator it = points.insert(make_pair(x, y)).first;

        if (bad(it)) {
            points.erase(it);
            return false;
        }

        area += point_area(it);

        while (it != points.begin() && bad(prev(it)))
            erase(prev(it));

        while (bad(next(it)))
            erase(next(it));

        return true;
    }

    // Returns 1 for strictly contained, 0 for on the border, and -1 for not contained.
    int contains(coord_t x, coord_t y) const {
        if (points.empty())
            return -1;

        map<coord_t, coord_t>::const_iterator first = points.begin(), last = prev(points.end());

        if (x < first->first || x > last->first)
            return -1;

        if (x == first->first)
            return y <= first->second ? 0 : -1;

        if (x == last->first)
            return y <= last->second ? 0 : -1;

        map<coord_t, coord_t>::const_iterator it = points.lower_bound(x);
        int64_t a = triangle_area(*prev(it), make_pair(x, y), *it);
        return a == 0 ? 0 : (a > 0 ? -1 : 1);
    }
};

// The combined hull with both upper and lower.
struct online_hull {
    upper_hull upper, lower;

    int64_t get_area_doubled() const {
        return upper.area + lower.area;
    }

    bool insert(coord_t x, coord_t y) {
        bool upper_insert = upper.insert(x, y);
        bool lower_insert = lower.insert(x, -y);
        return upper_insert || lower_insert;
    }

    int contains(coord_t x, coord_t y) const {
        int upper_contains = upper.contains(x, y);
        int lower_contains = lower.contains(x, -y);
        return min(upper_contains, lower_contains);
    }
};

online_hull hull;

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int type, x, y;

    while (cin >> type >> x >> y) {
        if (type == 1) {
            hull.insert(x, y);
            cout << hull.get_area_doubled() << '\n';
        } else if (type == 2) {
            int contains = hull.contains(x, y);
            cout << (contains == 0 ? "border" : (contains > 0 ? "inside" : "outside")) << '\n';
        } else {
            assert(false);
        }
    }

    cerr << "Hull size: " << hull.upper.size() + hull.lower.size() << '\n';
}
