#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

// Also known as "extended KMP"
template<typename T>
vector<int> z_algorithm(const T &pattern) {
    // Z[i] = for the suffix [i, n) of pattern, the longest prefix that is also a prefix of pattern.
    int n = int(pattern.size());
    vector<int> Z(n, 0);
    Z[0] = n;
    int loc = 1;

    for (int i = 1; i < n; i++) {
        if (i < loc + Z[loc])
            Z[i] = min(Z[i - loc], loc + Z[loc] - i);

        while (i + Z[i] < n && pattern[Z[i]] == pattern[i + Z[i]])
            Z[i]++;

        // Find the location with the furthest-reaching umbrella.
        if (i + Z[i] > loc + Z[loc])
            loc = i;
    }

    return Z;
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    string pattern, text;
    cin >> pattern >> text;
    int n = int(pattern.size()), m = int(text.size());

    vector<int> Z = z_algorithm(pattern + text);

    for (int i = 0; i < m; i++)
        if (Z[i + n] >= n)
            // Found a match starting at index i of text
            cout << i << '\n';
}
