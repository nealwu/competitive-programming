#include <algorithm>
#include <array>
#include <cassert>
#include <ctime>
#include <iomanip>
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
            dp[i + 1][j + 1] = max({dp[i + 1][j], dp[i][j + 1], dp[i][j] + (S[i] == T[j])});

    return dp[n][m];
}

int longest_common_subsequence(const string &S, const string &T) {
    int n = int(S.size());
    int m = int(T.size());
    vector<int> dp(m + 1, 0);

    for (int i = 0; i < n; i++) {
        vector<int> next_dp(m + 1, 0);

        for (int j = 0; j < m; j++)
            next_dp[j + 1] = max({next_dp[j], dp[j + 1], dp[j] + (S[i] == T[j])});

        dp.swap(next_dp);
    }

    return dp[m];
}

string construct_longest_common_subsequence_bitset(const string &S, const string &T) {
    int n = int(S.size());
    int m = int(T.size());
    vector<int> dp(m + 1, 0);
    vector<vector<bool>> previous(n + 1, vector<bool>(m + 1, false));

    for (int i = 0; i < n; i++) {
        vector<int> next_dp(m + 1, 0);

        for (int j = 0; j < m; j++) {
            next_dp[j + 1] = max({next_dp[j], dp[j + 1], dp[j] + (S[i] == T[j])});
            previous[i + 1][j + 1] = next_dp[j + 1] == next_dp[j];
        }

        dp.swap(next_dp);
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

string construct_longest_common_subsequence_hirschberg(const string_view &S, const string_view &T) {
    int n = int(S.size());
    int m = int(T.size());

    if (n == 0 || m == 0)
        return "";

    if (n == 1)
        return T.find(S[0]) == string::npos ? "" : string(1, S[0]);

    int mid = n / 2;
    vector<int> dp_first(m + 1, 0), dp_second(m + 1, 0);
    vector<int> next_dp(m + 1, 0);

    for (int i = 0; i < mid; i++) {
        for (int j = 0; j < m; j++)
            next_dp[j + 1] = max({next_dp[j], dp_first[j + 1], dp_first[j] + (S[i] == T[j])});

        dp_first.swap(next_dp);
    }

    next_dp.assign(m + 1, 0);

    for (int i = n - 1; i >= mid; i--) {
        for (int j = m - 1; j >= 0; j--)
            next_dp[j] = max({next_dp[j + 1], dp_second[j], dp_second[j + 1] + (S[i] == T[j])});

        dp_second.swap(next_dp);
    }

    int split = 0;

    for (int j = 1; j <= m; j++)
        if (dp_first[j] + dp_second[j] > dp_first[split] + dp_second[split])
            split = j;

    dp_first.clear();
    dp_second.clear();
    next_dp.clear();

    return (
        construct_longest_common_subsequence_hirschberg(S.substr(0, mid), T.substr(0, split)) +
        construct_longest_common_subsequence_hirschberg(S.substr(mid), T.substr(split))
    );
}

int main() {
    ios::sync_with_stdio(false);
#ifndef NEAL_DEBUG
    cin.tie(nullptr);
#endif

    cerr << setprecision(3);

    string S, T;
    cin >> S >> T;
    long double begin;

    begin = clock();
    int longest = longest_common_subsequence(S, T);
    cerr << "standard  : " << (clock() - begin) / CLOCKS_PER_SEC << 's' << endl;

    begin = clock();
    string answer = construct_longest_common_subsequence_bitset(S, T);
    cerr << "bitset    : " << (clock() - begin) / CLOCKS_PER_SEC << 's' << endl;

    begin = clock();
    assert(longest_common_subsequence_quadratic_memory(S, T) == longest);
    cerr << "quadratic : " << (clock() - begin) / CLOCKS_PER_SEC << 's' << endl;

    assert(int(answer.size()) == longest);
    assert(is_subsequence(answer, S) && is_subsequence(answer, T));

    begin = clock();
    string hirschberg_answer = construct_longest_common_subsequence_hirschberg(S, T);
    cerr << "hirschberg: " << (clock() - begin) / CLOCKS_PER_SEC << 's' << endl;

    assert(int(hirschberg_answer.size()) == longest);
    assert(is_subsequence(hirschberg_answer, S) && is_subsequence(hirschberg_answer, T));

    cout << longest << '\n';
    cout << answer << '\n';
}
