// https://en.wikipedia.org/wiki/Edit_distance
// https://en.wikipedia.org/wiki/Levenshtein_distance
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <vector>
using namespace std;

const int INF = int(1e9) + 5;

int edit_distance_quadratic_memory(const string &S, const string &T) {
    int n = int(S.size());
    int m = int(T.size());
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, INF));

    for (int i = 0; i <= n; i++)
        dp[i][0] = i;

    for (int j = 0; j <= m; j++)
        dp[0][j] = j;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            dp[i + 1][j + 1] = min({dp[i + 1][j] + 1, dp[i][j + 1] + 1, dp[i][j] + (S[i] != T[j])});

    return dp[n][m];
}

int edit_distance(const string &S, const string &T) {
    int n = int(S.size());
    int m = int(T.size());
    vector<int> dp(m + 1);
    iota(dp.begin(), dp.end(), 0);

    for (int i = 0; i < n; i++) {
        vector<int> ndp(m + 1, INF);
        ndp[0] = i + 1;

        for (int j = 0; j < m; j++)
            ndp[j + 1] = min({ndp[j] + 1, dp[j + 1] + 1, dp[j] + (S[i] != T[j])});

        dp.swap(ndp);
    }

    return dp[m];
}

vector<string> construct_edit_distance(const string &S, const string &T) {
    int n = int(S.size());
    int m = int(T.size());
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, INF));

    for (int i = 0; i <= n; i++)
        dp[i][0] = i;

    for (int j = 0; j <= m; j++)
        dp[0][j] = j;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            dp[i + 1][j + 1] = min({dp[i + 1][j] + 1, dp[i][j + 1] + 1, dp[i][j] + (S[i] != T[j])});

    vector<string> left = {S}, right = {T};

    while (n > 0 || m > 0) {
        if (n > 0 && dp[n][m] == dp[n - 1][m] + 1) {
            n--;
            string str = left.back();
            str.erase(str.begin() + n);
            left.push_back(str);
        } else if (m > 0 && dp[n][m] == dp[n][m - 1] + 1) {
            m--;
            string str = right.back();
            str.erase(str.begin() + m);
            right.push_back(str);
        } else if (n > 0 && m > 0 && dp[n][m] == dp[n - 1][m - 1] + (S[n - 1] != T[m - 1])) {
            n--;
            m--;

            if (S[n] != T[m]) {
                string str = left.back();
                str[n] = T[m];
                left.push_back(str);
            }
        } else {
            assert(false);
        }
    }

    assert(left.back() == right.back());
    right.pop_back();

    while (!right.empty()) {
        left.push_back(right.back());
        right.pop_back();
    }

    return left;
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    string S, T;
    cin >> S >> T;
    int dist = edit_distance(S, T);
    cout << dist << '\n';
    vector<string> answer = construct_edit_distance(S, T);

    for (string &str : answer)
        cout << str << '\n';

    assert(edit_distance_quadratic_memory(S, T) == dist);
    assert(int(answer.size()) == dist + 1);
    assert(answer.front() == S && answer.back() == T);

    for (int i = 0; i < dist; i++)
        assert(edit_distance(answer[i], answer[i + 1]) == 1);
}
