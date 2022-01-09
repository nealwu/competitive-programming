#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

// Use compare = less<T>() for a min heap and compare = greater<T>() for a max heap.
// When there are ties, the left value will be the parent of the right value.
template<typename T, typename Compare>
vector<int> build_cartesian_tree(const vector<T> &A, Compare &&compare) {
    int n = int(A.size());
    vector<int> parent(n, -1);
    vector<int> stack;

    for (int i = 0; i < n; i++) {
        int erased = -1;

        while (!stack.empty() && compare(A[i], A[stack.back()])) {
            erased = stack.back();
            stack.pop_back();
        }

        parent[i] = stack.empty() ? -1 : stack.back();

        if (erased >= 0)
            parent[erased] = i;

        stack.push_back(i);
    }

    return parent;
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    int N;
    cin >> N;
    vector<int> A(N);

    for (auto &a : A)
        cin >> a;

    vector<int> parent = build_cartesian_tree(A, less<int>());

    for (int i = 0; i < N; i++)
        cout << parent[i] + 1 << (i < N - 1 ? ' ' : '\n');

    parent = build_cartesian_tree(A, greater<int>());

    for (int i = 0; i < N; i++)
        cout << parent[i] + 1 << (i < N - 1 ? ' ' : '\n');
}
