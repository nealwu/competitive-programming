#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <vector>
using namespace std;

bool is_subsequence(const string &sub, const string &str) {
    size_t index = 0;

    for (char ch : str)
        if (index < sub.size() && ch == sub[index])
            index++;

    return index == sub.size();
}

int longest_common_subsequence_quadratic_memory(const string &S, const string &T) {
    int n = int(S.size());
    int m = int(T.size());
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            if (S[i] == T[j])
                dp[i + 1][j + 1] = dp[i][j] + 1;
            else
                dp[i + 1][j + 1] = max(dp[i][j + 1], dp[i + 1][j]);

    return dp[n][m];
}

int longest_common_subsequence(const string &S, const string &T) {
    int n = int(S.size());
    int m = int(T.size());
    vector<int> dp(m + 1, 0);

    for (int i = 0; i < n; i++) {
        vector<int> next_dp(m + 1, 0);

        for (int j = 0; j < m; j++)
            if (S[i] == T[j])
                next_dp[j + 1] = dp[j] + 1;
            else
                next_dp[j + 1] = max(dp[j + 1], next_dp[j]);

        swap(dp, next_dp);
    }

    return dp[m];
}

string construct_longest_common_subsequence(const string &S, const string &T) {
    int n = int(S.size());
    int m = int(T.size());
    vector<int> dp(m + 1, 0);
    vector<vector<bool>> previous(n + 1, vector<bool>(m + 1, false));

    for (int i = 0; i < n; i++) {
        vector<int> next_dp(m + 1, 0);

        for (int j = 0; j < m; j++)
            if (S[i] == T[j]) {
                next_dp[j + 1] = dp[j] + 1;
            } else {
                next_dp[j + 1] = max(dp[j + 1], next_dp[j]);
                previous[i + 1][j + 1] = next_dp[j + 1] == next_dp[j];
            }

        swap(dp, next_dp);
    }

    int a = n, b = m;
    string common;

    while (a > 0 && b > 0) {
        if (S[a - 1] == T[b - 1]) {
            common += S[a - 1];
            a--; b--;
            continue;
        }

        if (previous[a][b])
            b--;
        else
            a--;
    }

    reverse(common.begin(), common.end());
    return common;
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    string S, T;
    cin >> S >> T;
    int longest = longest_common_subsequence(S, T);
    string answer = construct_longest_common_subsequence(S, T);
    cout << longest << '\n';
    cout << answer << '\n';

    assert(longest_common_subsequence_quadratic_memory(S, T) == longest);
    assert(int(answer.size()) == longest);
    assert(is_subsequence(answer, S));
    assert(is_subsequence(answer, T));
}
