#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

namespace KMP {
    // Returns the next length to search from (the longest suffix of the current string that is a prefix of the pattern)
    // after starting from the prefix `len` and adding char `c`.
    // Runs in worst case O(len) but amortizes to O(1) in most situations.
    template<typename T, typename T_elem>
    int get_link(const T &pattern, const vector<int> &fail, int len, const T_elem &c) {
        while (len > 0 && pattern[len] != c)
            len = fail[len];

        if (pattern[len] == c)
            len++;

        return len;
    }

    // Computes the failure function on `pattern` so that we can search for it in future strings.
    template<typename T>
    vector<int> compute_failure_function(const T &pattern) {
        // fail[i] = for the prefix [0, i) of `pattern`, the length of the longest proper prefix that is also a suffix.
        int n = int(pattern.size());
        vector<int> fail(n + 1, 0);
        int len = 0;

        for (int i = 1; i < n; i++) {
            len = get_link(pattern, fail, len, pattern[i]);
            fail[i + 1] = len;
        }

        return fail;
    }

    template<typename T>
    vector<int> find_matches(const T &pattern, const T &text, const vector<int> &fail) {
        if (pattern.size() > text.size())
            return {};

        vector<int> matches;
        int n = int(pattern.size()), m = int(text.size());
        int len = 0;

        for (int i = 0; i < m; i++) {
            len = get_link(pattern, fail, len, text[i]);

            // Check for a match whose last character is at index i.
            if (len == n) {
                matches.push_back(i - (n - 1));
                len = fail[len];
            }
        }

        return matches;
    }

    // Finds all indices where `pattern` occurs within `text`.
    template<typename T>
    vector<int> find_matches(const T &pattern, const T &text) {
        return find_matches(pattern, text, compute_failure_function(pattern));
    }
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    string pattern, text;
    cin >> pattern >> text;
    vector<int> matches = KMP::find_matches(pattern, text);

    for (int match : matches)
        cout << match << '\n';
}
