#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

const int INF = int(1e9) + 5;

struct point {
    int64_t x, y;
    int index;

    point() : x(0), y(0), index(-INF) {}

    point(int64_t _x, int64_t _y, int _index = -INF) : x(_x), y(_y), index(_index) {}

    point& operator+=(const point &other) { x += other.x; y += other.y; return *this; }
    point& operator-=(const point &other) { x -= other.x; y -= other.y; return *this; }

    point operator+(const point &other) const { return point(*this) += other; }
    point operator-(const point &other) const { return point(*this) -= other; }

    bool operator==(const point &other) const { return x == other.x && y == other.y; }
    bool operator!=(const point &other) const { return !(*this == other); }

    point operator-() const {
        return point(-x, -y);
    }

    bool top_half() const {
        return y > 0 || (y == 0 && x > 0);
    }

    int64_t norm() const {
        return (int64_t) x * x + (int64_t) y * y;
    }

    double dist() const {
        return sqrt(norm());
    }

    friend ostream& operator<<(ostream &os, const point &p) {
        return os << '(' << p.x << ", " << p.y << ')';
    }
};

int64_t manhattan_dist(const point &a, const point &b) {
    return (int64_t) abs(a.x - b.x) + abs(a.y - b.y);
}

// Sort in increasing order of y, with ties broken in increasing order of x.
bool yx_compare(const point &a, const point &b) {
    return make_pair(a.y, a.x) < make_pair(b.y, b.x);
}

struct union_find {
    // When data[x] < 0, x is a root and -data[x] is its tree size. When data[x] >= 0, data[x] is x's parent.
    vector<int> data;
    int components = 0;

    union_find(int n = -1) {
        if (n >= 0)
            init(n);
    }

    void init(int n) {
        data.assign(n + 1, -1);
        components = n;
    }

    int find(int x) {
        return data[x] < 0 ? x : data[x] = find(data[x]);
    }

    int get_size(int x) {
        return -data[find(x)];
    }

    bool unite(int x, int y) {
        x = find(x);
        y = find(y);

        if (x == y)
            return false;

        if (-data[x] < -data[y])
            swap(x, y);

        data[x] += data[y];
        data[y] = x;
        components--;
        return true;
    }
};


struct edge {
    int index1, index2;
    int64_t dist;

    edge() {}

    edge(int _index1, int _index2, int64_t _dist) : index1(_index1), index2(_index2), dist(_dist) {}

    bool operator<(const edge &other) const {
        return dist < other.dist;
    }
};

point rotate90(point p) {
    swap(p.x, p.y);
    p.x = -p.x;
    return p;
}

void rotate_all(vector<point> &points) {
    for (point &p : points)
        p = rotate90(p);
}

void swap_all(vector<point> &points) {
    for (point &p : points)
        swap(p.x, p.y);
}

bool has_better_sum(const point &a, const point &b) {
    if (a.index < 0)
        return false;

    if (b.index < 0)
        return true;

    return a.x + a.y < b.x + b.y;
}

vector<point> buffer, c_buffer;

void solve(vector<point> &points, vector<point> &closest, int start, int end) {
    if (end - start <= 1)
        return;

    int mid = (start + end) / 2;
    solve(points, closest, start, mid);
    solve(points, closest, mid, end);
    int right = mid, n = 0;
    point min_sum;

    // Merge sort by y - x, and keep track of the smallest x + y point we've seen so far.
    // Thus for each left point, we find the right point with the minimum value of x + y that satisfies
    // yx_compare(left, right) = true and right.y - right.x <= left.y - left.x.
    for (int i = start; i < mid; i++) {
        while (right < end && points[right].y - points[right].x <= points[i].y - points[i].x) {
            buffer[n] = points[right];
            c_buffer[n] = closest[right];

            if (has_better_sum(points[right], min_sum))
                min_sum = points[right];

            n++;
            right++;
        }

        if (has_better_sum(min_sum, closest[i]))
            closest[i] = min_sum;

        buffer[n] = points[i];
        c_buffer[n] = closest[i];
        n++;
    }

    // Copy the results back in their sorted order.
    copy(buffer.begin(), buffer.begin() + n, points.begin() + start);
    copy(c_buffer.begin(), c_buffer.begin() + n, closest.begin() + start);
}

vector<edge> manhattan_mst(vector<point> points) {
    int N = int(points.size());
    vector<point> closest;
    vector<edge> edges;

    buffer.resize(N);
    c_buffer.resize(N);

    // We find four edges for each point: one to the closest point in each of the four octants on the point's right.
    for (int rep = 0; rep < 4; rep++) {
        sort(points.begin(), points.end(), yx_compare);
        closest.assign(N, point());

        // For each point (x, y), find the closest point (x', y') such that y' >= y and y' - x' <= y - x.
        // In other words, this finds the closest point in the octant to the right of the point and slightly above.
        // See https://www.topcoder.com/community/competitive-programming/tutorials/line-sweep-algorithms/ for details.
        solve(points, closest, 0, N);

        for (int i = 0; i < N; i++)
            if (closest[i].index >= 0)
                edges.emplace_back(points[i].index, closest[i].index, manhattan_dist(points[i], closest[i]));

        if (rep % 2 == 0)
            swap_all(points);
        else
            rotate_all(points);
    }

    // To return the points back to normal:
    // rotate_all(points); rotate_all(points);

    sort(edges.begin(), edges.end());
    union_find UF(N);
    vector<edge> mst;

    for (edge &e : edges)
        if (UF.unite(e.index1, e.index2))
            mst.push_back(e);

    return mst;
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N;
    cin >> N;
    vector<point> points(N);

    for (point &p : points)
        cin >> p.x >> p.y;

    for (int i = 0; i < N; i++)
        points[i].index = i;

    vector<edge> mst = manhattan_mst(points);
    assert(int(mst.size()) == max(N - 1, 0));
    int64_t total = 0;

    for (edge &e : mst)
        total += e.dist;

    cout << total << '\n';
}
